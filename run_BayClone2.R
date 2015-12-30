run_BayClone2 <- function(random_p0,with_M)
{
    #### with_M can be 1 or 0;  random_p0 can be 1 or 0
    #### when T = 1, default setting is random_p0 = 0    
    set.seed(79861)
    #library(BayClone2)
    library(combinat)
    source("./plotResults.R")
    
    source("./fn_BayClone2_all.R")
    dyn.load("./BayClone2_all.so")
    

    ######################################
    S <- 100 # THE NUMBER OF SNP (I.E. NUMBER OF ROWS OF Z)
    T <- 4 #THE NUMBER OF TISSUES
    
    
    ########################################################
    #SIMULATE DATA
    ########################################################
    #DATA
    M <- array(NA, dim=c(S, T)) #SAMPLE COPY NUMBER FROM BIO-X CALLER
    N <- array(NA, dim=c(S, T)) #NUMBER OF READS INDICATING MUTATIONS
    n <- array(NA, dim=c(S, T)) #TOTAL NUMBER OF READS
    
    
    #SIMULTATION TRUTH
    sim.tr <- NULL
    sim.tr$C <- 3 #INCLUDING BACKGROUND
    sim.tr$const <- 5
    
    
    #SET L AND Z
    sim.tr$L <- array(2, dim=c(S, sim.tr$C))
    sim.tr$L[45:50,2] <- 3
    sim.tr$L[51:55,2] <- 1
    
    sim.tr$L[,3] <- sim.tr$L[,2]
    sim.tr$L[71:75,3] <- 0
    sim.tr$L[91:95,3] <- 1
    sim.tr$L[96:100,3] <- 0
    
    
    
    sim.tr$Z <- array(0, dim=c(S, sim.tr$C))
    sim.tr$Z[,1] <- 2
    
    sim.tr$Z[56:60,2] <- 1
    sim.tr$Z[66:80,2] <- 1
    sim.tr$Z[81:90,2] <- 2
    
    
    sim.tr$Z[,3] <- sim.tr$Z[,2]
    sim.tr$Z[16:20, 3] <- 1
    sim.tr$Z[21:30, 3] <- 2
    sim.tr$Z[31:40, 3] <- 1
    sim.tr$Z[71:75,3] <- 0
    sim.tr$Z[96:100, 3] <- 0
    
    
    cbind(sim.tr$L, sim.tr$Z)
    sum(sim.tr$Z > sim.tr$L)
    
    ##PHI
    sim.tr$phi <- rgamma(T, 600, 3)  #TOTAL NUMBER OF READS
    
    
    ###W--WEIHTS
    sim.tr$w <- sim.tr$w_star <- array(NA, dim=c(T, sim.tr$C))
    a.vec <- c(0.01, 0.75, 0.25)*40
    
    for(i.t in 1:T)
    {
        sim.tr$w_star[i.t, ] <- rgamma(sim.tr$C, a.vec, 1)
        sim.tr$w[i.t,] <- sim.tr$w_star[i.t,]/sum(sim.tr$w_star[i.t,])
    }
    
    sim.tr$mu <- sim.tr$p <- array(NA, dim=c(S, T))
    
    sim.tr$p_z0 <- c(0.05, rep(1, sim.tr$C-1))
    
    for(i.s in 1:S)
        for(i.t in 1:T)
        {
            sim.tr$mu[i.s, i.t] <- sum(sim.tr$w[i.t,]*sim.tr$L[i.s,])
            sim.tr$p[i.s, i.t] <- sum(sim.tr$w[i.t,]*sim.tr$Z[i.s,]*sim.tr$p_z0)/sum(sim.tr$w[i.t,]*sim.tr$L[i.s,])
        }
    
    
    N <- n <- array(0, dim=c(S, T))
    for(i.t in 1:T)
    {
        N[,i.t] <- rpois(S, sim.tr$phi[i.t]*sim.tr$mu[,i.t]/2)  ##LET'S LEAVE THE GENERATION OF N_ST AS IT IS NOW (NOTHING WILL BE HURTED.
        n[,i.t] <- rbinom(S, N[,i.t], sim.tr$p[,i.t])
    }
    
    M <- log(sim.tr$mu/2, 2)
    
    #library("BayClone2")
    
    ##READ IN DATA
    #data(BayClone2_Simulation1_tot)
    #data(BayClone2_Simulation1_mut)
    #data(BayClone2_Simulation1_log_M)
    ##TOTAL NUMBER OF READS AT LOCUS s IN SAMPLE t
    #N <- as.matrix(BayClone2_Simulation1_tot)  
    ##NUMBER OF READS WITH VARIANT SEQUENCE AT LOCUS s IN SAMPLE t
    #n <- as.matrix(BayClone2_Simulation1_mut) 
    ### log2 ratio of CNV 
    #M <- as.matrix(BayClone2_Simulation1_log_M)
    
    S <- nrow(N)  # THE NUMBER OF LOCI (I.E. NUMBER OF ROWS OF N (AND n))
    T <- ncol(N) #THE 
    
    print(paste("M has dim ", dim(M)))
    print(paste("N has dim ", dim(N)))
    print(paste("n has dim ", dim(n)))
    
    ###################################
    ## BAYCLONE2 *WITHOUT* M AND WITH RANDOM P0
    ###################################
    #####------ DATA
    S #--NUMBER OF LOCI (NUMBER OF ROWS IN THE DATA MATRIX)
    T #--NUMBER OF SAMPLES (NUMBER OF ROWS IN THE DATA MATRIX)
    N ##--TOTAL COUNTS
    n ##--COUNTS WITH VARIANTS
    
    #HYPER-PARAMETER
    hyper <- NULL
    
    #NUMBER OF CELL TYPES (GEOMETRIC DIST)
    hyper$r <- 0.2
    
    #PRIOR FOR cIBP
    hyper$Q <- 3  #NUMBER OF COPIES -- q = 0, 1, 2, 3
    ##BETA-DIRICHLET
    hyper$alpha <- 2
    hyper$beta <- 1
    hyper$gam <- rep(0.5, hyper$Q)
    
    #PRIOR FOR PO
    hyper$a_z0 <- 0.3
    hyper$b_z0 <- 5
    hyper$p0 <- 0.05
    
    #PRIOR FOR PHI--TOTAL NUMBER OF READS IN A SAMPLE T
    hyper$b <- 3
    hyper$a <- median(N)*hyper$b
    
    #PRIOR FOR W
    hyper$d0 <- 0.5
    hyper$d <- 1
    
    ###PRIOR FOR M
    hyper$sig2 <- rep(0.00125, T)
    
    #########----MCMC PARAMETERS
    n.sam <- 500; burn.in <- 100
    Min_C <- 2  ##INCLUDING THE BACKGROUND SUBCLONE (SUBCLONE 0)
    Max_C <- 16   ##INCLUDING THE BACKGROUND SUBCLONE (SUBCLONE 0)
    ave.B <- 0.025 ###--a portion for training data set
    
    #random_p0 = 0;
    #with_M = 0

    MCMC.sam <- BayClone2(Min_C, Max_C, S, T, burn.in, n.sam, N, n, hyper, ave.B, random_p0, with_M, MM=M)

    dir.create("OUT")
    dir.create("PLOTS")
    #### save those samples
    if(is.numeric(MCMC.sam) && (MCMC.sam == -1))
    {
        print("did not run !!!")
    }
    else
    {
        rdaMCMCFile = sprintf("./OUT/MCMC_random_p0_%d_with_M%d.rda",random_p0,with_M);
        print(paste("MCMC filename ",rdaMCMCFile))
        save(MCMC.sam,file=rdaMCMCFile)
    
        print("all done !!")    
        
        xt = tabulate(MCMC.sam$C)
        postMode_C = which(xt==max(xt))
        
         
        post_dist_C <- fn_post_C(MCMC.sam$C, Min_C, Max_C)
        
        CC <- postMode_C
        
        print(paste0("Posterior Mode of C is ",postMode_C))
        
        point.est <- fn_posterior_point(CC, S, T, MCMC.sam,random_p0,with_M)
        
        rdaPostFile = sprintf("./OUT/Post_random_p0_%d_with_M%d_%d.rda",random_p0,with_M,CC)
        save(post_dist_C,point.est,file=rdaPostFile)
        
        print("rda files are generated !!")
        print("NOW PLOTTING .... ")
    
        plotResults(random_p0,with_M,N,n)
        print("Plotting done !!!")    
    
    }
    
}

args <- commandArgs(TRUE)
vcf_fileName = args[1]
cnv_fileName = args[2]

run_BayClone2(1,1)

print(paste0("vcf file: ",vcf_fileName," cnv file: ",cnv_fileName))
