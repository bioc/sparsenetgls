#'The glassonet2() function
#'@description The glassonet2 function is designed to learn the graph 
#'structure, the corresponding precision matrix and covariance matrix 
#'by using the graph lasso method.
#'
#'@import huge Matrix MASS
#'
#'@param Y0 The data matrix for the GGM model.
#'
#'@param nlambda The number of interval used in the penalized path in lasso
#'and elastics. It results in the number of lambda values to be used in the
#'penalization. The default value is nlambda assigned in the parent function
#'sparsenetgls().
#'
#'@param lambda.min.ratio It is the default parameter set in function huge()
#'in the package 'huge'. Quoted from huge(), it is the minimal value of lambda,
#'being a fraction of the upper bound (MAX) of the regularization/ thresholding
#'parameter that makes all the estimates equal to 0.
#'The default value is 0.001.
#'
#'@param method There are two options for the method parameter which is
#'provided in the huge() function. One is 'glasso' and the other one is 'mb'.
#'
#'@return Return the precision matrix 'OMEGAMATRIX', penalized path parameter
#'lambda 'lambda' and covariance matrix 'COVMATRIX'.  
#'
#'
#'@examples
#'n=20
#'VARknown <- rWishart(1, df=4, Sigma=matrix(c(1,0,0,0,1,0,0,0,1),
#'nrow=3,ncol=3))
#'Y0 <- mvrnorm(n=n,mu=rep(0.5,3),Sigma=VARknown[,,1])
#'fitglasso <- glassonet2(Y0=Y0,nlambda=5,method='glasso')
#'
#'@export

glassonet2 <- function(Y0, nlambda = nlambda, lambda.min.ratio = 0.001, 
    method) {
    p <- dim(Y0)[2]
    n <- dim(Y0)[1]
    solg <- huge(Y0, method = method, nlambda = nlambda, 
        lambda.min.ratio = lambda.min.ratio, scr = TRUE, 
        cov.output = TRUE)
    
    COVMATRIX = OMEGAMATRIX <- array(dim = c(p, p, 
        nlambda))
    
    for (i in seq_len(nlambda)) {
        OMEGAMATRIX[, , i] <- as.matrix(solg$icov[[i]])
        COVMATRIX[, , i] <- as.matrix(solg$cov[[i]])
    }
    
    return(list(OMEGAMATRIX = OMEGAMATRIX, lambda = solg$lambda, 
        COVMATRIX = COVMATRIX))
}
