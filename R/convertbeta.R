#'The convertbeta() function
#'@description The covertbeta function is designed to convert the regression
#'coefficients derived from the standardized data.
#'
#'@param X It is a dataset of explanatory variables.
#'
#'@param Y It is the multivariate response variables. 
#'
#'@param q It is an integer representing the number of explanatory variables.
#'
#'@param beta0 The vector contains the regression coefficients from
#'sparsenetgls. 
#'
#'@return Return the list of converted regression coefficients of the
#'explanatory variables 'betaconv' and intercept value 'betaconv_int'. 
#'
#'
#'@examples
#'X <- mvrnorm(n=20,mu=rep(0,5),Sigma=Diagonal(5,rep(1,5)))
#'Y <- mvrnorm(n=20,mu=rep(0.5,10),Sigma=Diagonal(10,rep(1,10)))
#'fitmodel <-  sparsenetgls(responsedata=Y,predictdata=X,nlambda=5,ndist=2,
#'method='elastic')
#'convertbeta(X=X,Y=Y,q=5+1,beta0=fitmodel$beta[,1])
#'
#'@export

convertbeta <- function(X, Y, q, beta0) {
    betaconv = beta0
    meanX <- colMeans(X, na.rm = TRUE)
    meanY <- mean(Y, na.rm = TRUE)
    sdX <- sqrt(rowVar(t(X)))
    sdY <- apply(Y, 2, sd, na.rm = TRUE)
    conv_ratio <- (mean(sdY) * sdX^(-1))
    betaconv_int <- -conv_ratio * meanX * beta0[2:q] + 
        beta0[1] * mean(sdY) + meanY
    betaconv[1] <- mean(betaconv_int, na.rm = TRUE)
    betaconv[2:q] <- conv_ratio * beta0[2:q]
    return(list(betaconv = betaconv, betaconv_int = betaconv_int))
}
