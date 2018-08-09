## ----setup, include=TRUE---------------------------------------------------
library(knitr)
library(rmarkdown)
knitr::opts_chunk$set(echo = TRUE)

## ----installation1,eval=FALSE----------------------------------------------
#  source("https://bioconductor.org/biocLite.R")
#  BiocInstaller::biocLite("sparsenetgls")
#  library(sparsenetgls)

## ----installation2,eval=FALSE----------------------------------------------
#  devtools::install_github("superOmics/sparsenetgls")
#  library(sparsenetgls)

## ----example 1: assess estimated result from a GGM, eval=TRUE--------------
library(sparsenetgls)
library(MASS)
library(Matrix)
#simulate the dataset
data(bandprec, package="sparsenetgls")
varKnown <- solve(as.matrix(bandprec))
prec <- as.matrix(bandprec)
    
Y0 <- mvrnorm(n=100, mu=rep(0,50), Sigma=varKnown)
nlambda=10                                     
#u-beta
u <- rep(1,8)
X_1 <- mvrnorm(n=100, mu=rep(0,8), Sigma=Diagonal(8,rep(1,8)))
Y_1 <- Y0+as.vector(X_1%*%as.matrix(u))
databand_X=X_1
databand_Y=Y_1

#produce the precision matrices
omega <- sparsenetgls(responsedata=databand_Y, predictdata=databand_X, nlambda=10,
ndist=1, method="glasso")$PREC_seq
omega_est <- array(dim=c(50,50,10))
for (i in seq_len(10)) omega_est[,,i] <- as.matrix(omega[[i]])

roc_path_result <- path_result_for_roc(PREC_for_graph=prec, OMEGA_path=omega_est
, pathnumber=10)
plot_roc(result_assessment=roc_path_result, group=FALSE, ngroup=0,
est_names="glasso estimation")

## ----example 2: Use the sparsenetgls and convertbeta function,eval=TRUE----

fitgls <- sparsenetgls(databand_Y, databand_X, nlambda=10, ndist=5, method="glasso")

#Convert the regression coefficients to its original scale
q <- dim(databand_X)[2]
nlambda=10

betagls <- matrix(nrow=nlambda, ncol=q+1)
for (i in seq_len(nlambda))   
betagls[i,] <- convertbeta(Y=databand_Y, X=databand_X, q=q+1,
beta0=fitgls$beta[,i])$betaconv

#Beta selection 

#select lamda and dist value based on the minimal variance of beta
ndist <- max(fitgls$power)-1
tr_gamma <- matrix(nrow=10, ncol=ndist-1)

for (j in seq_len(ndist-1))
    for (i in seq_len(nlambda)) 
        tr_gamma[i,j] <- (sum(diag(fitgls$covBeta[,,j,i]))) 

select.lambda.dist <- which(tr_gamma==min(tr_gamma), arr.ind=TRUE)
select.lambda.dist
betagls_select <- betagls[select.lambda.dist[1],]

#row is lambda and column is dist
varbeta <- diag(fitgls$covBeta[,,ndist,select.lambda.dist[1]])

##select lamda and dist value based on the AIC and BIC
select.lambda.dist2 <- which(fitgls$bic==min(fitgls$bic,na.rm=TRUE), arr.ind=TRUE)
select.lambda.dist3 <- which(fitgls$aic==min(fitgls$aic,na.rm=TRUE), arr.ind=TRUE)
varbeta_bic <- diag(fitgls$covBeta[,,ndist,select.lambda.dist2[1]])
varbeta_aic <- diag(fitgls$covBeta[,,ndist,select.lambda.dist3[1]])


## ----example 3, Get the visualized result, eval=TRUE-----------------------
plotsngls(fitgls,ith_lambda=5)

## ----example 4, Use different options of GGM estimation, eval=TRUE---------
#Use the glasso method to estimate the precision matrix
fitgls_g <- sparsenetgls(databand_Y, databand_X, nlambda=10, ndist=5, method="elastic")

#Uset the lasso method to approximate the precision matrix
#fitgls_l <- sparsenetgls(databand_Y, databand_X, nlambda=10, ndist=5, method="lasso")

#use the Meinshausen B?hlmann method to approximate the precision matrix
#fitgls_m <- sparsenetgls(databand_Y, databand_X, nlambda=10, ndist=5, method="mb")

## ----sessionInfo,echo=TRUE,eval=TRUE---------------------------------------
sessionInfo()

