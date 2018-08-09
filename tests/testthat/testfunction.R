#'Testthat program for packaging
#'

test_that("The sparsenetgls() using a diagonal precision matrix",
{
 set.seed(100)
 theta <- rgamma(10,1,1)^2
 varknown <- Diagonal(10,theta)
 prec <- solve(varknown)
 y <- mvrnorm(500,mu=rep(0,10),Sigma=varknown)
 beta <- c(0,1,2,3,4)
 x <- mvrnorm(500,mu=rep(0,5),Sigma=Diagonal(5,rep(1,5)))
 y1 <- as.vector(x%*%beta)+y

 fitsparse <- sparsenetgls(responsedata=y1,lambda.min.ratio=1e-7,predictdata=x,nlambda=5,ndist=1,method="glasso")

 #conversion to the original scale
 beta_gls <- convertbeta(Y=y1,X=x,q=6,beta0=fitsparse$beta[,1])$betaconv

 coeff=matrix(nrow=10,ncol=5)
 for ( i in 1:10) {lmf=lm(y1[,i]~x)
 coeff[i,]=lmf$coefficients[2:6]
 }
 lmbeta0=colMeans(coeff)

 expect_equal(beta_gls[2:6],lmbeta0,tolerance = 0.001,scale=1)
})

test_that("test the dist_tune function",{
          set.seed(200)
          y <- mvrnorm(n=100,mu=rep(0,5),Sigma=Matrix::Diagonal(5,rep(1,5)))
          covy=cov(y)
          covy
          covy[2,3]=covy[3,2]=covy[2,5]=covy[5,2]=covy[3,4]=covy[4,3]=covy[4,5]=covy[5,4]=0
          prec=solve(covy)
          omega_seq=huge::huge(covy,lambda.min.ratio=0.00001,nlambda=10,scr=TRUE,method="glasso")
          cov_sele=dist_tune(covy=covy,covstart=omega_seq$icov[[5]],ndist=7,p=5)
          powers=cov_sele$power
          expect_equal(cov_sele$cov_adj[,,powers],covy,scale=1)
}
)










