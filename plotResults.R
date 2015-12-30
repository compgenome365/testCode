#################################################################        
plotResults <- function(random_p0,with_M,N,n)
{    
    library(gplots);
    S <- nrow(N)  # THE NUMBER OF LOCI (I.E. NUMBER OF ROWS OF N (AND n))
    T <- ncol(N)  #THE NUMBER OF TISSUE SAMPLES  (I.E. NUMBER OF COLUMNS OF N (AND n))
    
    #sampleID = 1
    Min_C <- 2  ##INCLUDING THE BACKGROUND SUBCLONE
    Max_C <- 16  ##INCLUDING THE BACKGROUND SUBCLONE
    
    
    mcmcFileName = sprintf("./OUT/MCMC_random_p0_%d_with_M%d.rda",random_p0,with_M);
    load(mcmcFileName)  
    
    xt = tabulate(MCMC.sam$C)
    postMode_C = which(xt==max(xt))
    CC <- postMode_C
    
    rdaPostFile = sprintf("./OUT/Post_random_p0_%d_with_M%d_%d.rda",random_p0,with_M,CC)
    load(rdaPostFile)  
    
    
    
    
    ####SUMMARY FOR FIGURES######################
    ##PLOT THE MARGINAL DISTRIBUTION OF C
    print("plotting C")
    pdfFileName = sprintf("./PLOTS/post_C_random_p0_%d_with_M%d_%d.pdf",random_p0,with_M,CC)
    pdf(pdfFileName)
    par(mar=c(4.5, 4.5, 2.1, 2.1))
    plot(post_dist_C[,1], post_dist_C[,2], type="o", lwd=4, col=1,xlab="C", ylab="P(C|n, N)", main="", cex.axis=1.5, cex.lab=1.5)
    abline(v=CC-1, lty=2, lwd=3, col=2)  #ASSUMED TRUE C=2
    dev.off()
    
    if(CC > 2)
    {
        #################################
        ### posterior point estimate of L
        print("plotting L")
        pdfFileName = sprintf("./PLOTS/post_L_random_p0_%d_with_M%d_%d.pdf",random_p0,with_M,CC)
        pdf(pdfFileName)
        lmat = rbind(c(0,3),c(2,1),c(0,4))
        lwid = c(1.5,4)
        lhei = c(1.5,4,1.2)
        post.L <- point.est$L
        colnames(post.L) <- (0:(CC-1))
        heatmap.2(post.L[,-1], trace="none", Colv=FALSE, Rowv=FALSE,dendrogram="none", scale="none", col=redgreen, colsep=(1:(CC-1)),sepcol=c("white", "white"), sepwidth=c(0.05, 0.1),main="",lmat = lmat, lwid = lwid, lhei = lhei)
        dev.off()
        
        #################################
        ### posterior point estimate of Z
        print("plotting Z")
        pdfFileName = sprintf("./PLOTS/post_Z_random_p0_%d_with_M%d_%d.pdf",random_p0,with_M,CC)
        pdf(pdfFileName)
        post.Z <- point.est$Z
        colnames(post.Z) <- (0:(CC-1))
        heatmap.2(post.Z[,-1], trace="none", Colv=FALSE, Rowv=FALSE,dendrogram="none", scale="none", col=redgreen, colsep=(1:(CC-1)),sepcol=c("white", "white"), sepwidth=c(0.05, 0.1),main="",lmat = lmat, lwid = lwid, lhei = lhei)
        dev.off()
    }
	else
	{
		printf("No plot for L and Z !!");
	}
    #################################
    ### posterior point estimate of w
    print("Writing C")
    post.w <- point.est$w
    txtFileName = sprintf("./PLOTS/post_w_random_p0_%d_with_M%d_%d.txt",random_p0,with_M,CC)
    write.table(post.w,txtFileName,sep="\t")
    if(CC > 2 && T > 1)
    {    
        print("plotting w")
        pdfFileName = sprintf("./PLOTS/post_w_random_p0_%d_with_M%d_%d.pdf",random_p0,with_M,CC)
        pdf(pdfFileName)
        colnames(post.w) <- (0:(CC-1))
        heatmap.2(post.w, trace="none", Colv=FALSE, Rowv=TRUE, dendrogram="none", scale="none", col=redgreen, colsep=(1:(CC-1)), sepcol=c("white", "white"), sepwidth=c(0.05, 0.1), main="",lmat = lmat, lwid = lwid, lhei = lhei)
        dev.off()
    }
	else
	{
		printf("No Plot for w !!");
	}
    ##################################
    nBreak = 50; #### for histogram
    
    ### PLOT N and n
    print("plotting N")
    pdfFileName = sprintf("./PLOTS/N_tot_random_p0_%d_with_M%d_%d.pdf",random_p0,with_M,CC)
    pdf(pdfFileName)
    hist(N,breaks=nBreak,col='blue1');
    dev.off();
    
    print("plotting n/N")
    pdfFileName = sprintf("./PLOTS/n_N_raio_random_p0_%d_with_M%d_%d.pdf",random_p0,with_M,CC)
    pdf(pdfFileName)
    hist(n/N,breaks=nBreak,xlim=c(-1.0,1.0),col='blue1');
    dev.off();
    
    #################################
    #### model checking plot
    ratio_1 = n/N;
    p_true = ratio_1[,1]
    p_est = MCMC.sam$p[,,1];
    p_est_mean = colMeans(p_est);
    
    diff_p = p_true-p_est_mean;
    
    #attach(mtcars)
    #par(mfrow=c(2,1)) 
    print("plotting diff_p")
    pdfFileName = sprintf("./PLOTS/p_diff_random_p0_%d_with_M%d_%d.pdf",random_p0,with_M,CC)
    pdf(pdfFileName)
    hist(diff_p,breaks=nBreak,xlim=c(-0.5,0.5),col='blue1');
    dev.off();
    
}
