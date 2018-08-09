#'The plot_roc() function is designed to produce the Reciever Operative
#'Characteristics (ROC) Curve for visualizing
#'the prediction accuracy of a Gaussian Graphical model (GGM) to the true
#'graph structure. The GGM must use a l-p
#'norm regularizations (p=1,2) with the series of solutions conditional on the
#'regularization parameter.
#'
#'@import Matrix MASS
#'@param result_assessment It is the list result from function
#'path_result_for_roc() which has five-dimensions recording
#'the path number (i.e. the order of \eqn{\lambda} ), the sensitivity,
#'the specificity, the Negative predicted value (NPV) and
#'the Postive predicted value (PPV) respectively.
#'
#'@param group It is a logical parameter indicating if the result_assessment
#'is for several GGM models. When it is TRUE,
#'it produceS the ROC from several GGM models. when it is FALSE, it only
#'produces a ROC for one model.
#'
#'@param ngroup It is an interger recording the number of models when group
#'is TRUE.
#'
#'@param est_names it is used for labeling the GGM model in legend of ROC
#'curve.
#'
#'@return  Return the plot of Receiver Operational Curve
#'
#'@examples
#'prec1 <- matrix(c(0,2,3,1,0,0.5,0,0,0.4),nrow=3,ncol=3)
#'Omega_est <- array(dim=c(3,3,3))
#'Omega_est[,,1] <- matrix(c(1,1,1,0.2,0.5,0.2,2,0.2,0.3),nrow=3,ncol=3)
#'Omega_est[,,2] <- matrix(c(0,1,1,1,0,0,0,0,1),nrow=3,ncol=3)
#'Omega_est[,,3] <- matrix(c(0,0,0,0,0,0,0,0,0),nrow=3,ncol=3)
#'roc_path_result <- path_result_for_roc(PREC_for_graph=prec1,
#'OMEGA_path=Omega_est,pathnumber=3)
#'plot_roc(result_assessment=roc_path_result,group=FALSE,ngroup=0,
#'est_names="test example")
#'
#'
#'@export

plot_roc <- function(result_assessment,group=TRUE,ngroup=0,est_names){

if (!group)
{Sensitivity <- result_assessment[,1]
FPR <- 1-result_assessment[,2]

plot(FPR,Sensitivity,xlab="1-specificity",ylab="sensitivity",ylim=c(0,1),
    xlim=c(0,1),pch=19,cex=1.5,main="ROC of network structure estimates")
lines(FPR,Sensitivity)

}  else {
    
    k=1;sentivity=FPR=vector("numeric",length=ngroup)
    
    plot((1-result_assessment[,2,1]),result_assessment[,1,1],
    xlab="1-specificity",ylab="sensitivity",
    pch=19,cex=1.5,ylim=c(0,1),xlim=c(0,1),
    main="ROC of network structure estimates")
    
    while (k<=ngroup){
    sensitivity=result_assessment[,1,k]
    FPR=1-result_assessment[,2,k]
    points(FPR,sensitivity,xlab="1-specificity",ylab="sensitivity",pch=19,
    cex=1,col=k)
    lines(FPR,sensitivity,col=k)
    k=k+1
    }
}
    segments(0,0,1,1)
    legend(0.1,0.90,est_names,lty=seq_len(ngroup),col=seq_len(ngroup),bty="n")
}
