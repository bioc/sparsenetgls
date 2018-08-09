#'The lassoglmnet() function is designed to learn the graph structure by
#'using the lasso and elastics net method.
#'@import glmnet Matrix MASS
#'
#'@param Y0 The data matrix for the GGM model.
#'
#'@param nlambda The number of interval used in the penalized path in lasso
#'and elastics. It results in the number of lambda values to be used in the
#'penalization. The default value is 10.
#'
#'@param alpha The vaule to be used in enet, it has values betwee 0 and 1.
#'The value of 0 is corresponding to l-1 penalization,
#'and 1 is corresponding to the l-2 regularization (Ridge regression).
#'The other values between 0 and 1 will result in a
#'combination of l1-l2 norm regularization named as elastic net.
#'
#'@return Return the regression coefficients of glmnet "coef_glmnet",
#'residuals from the glmnet "resid_glmnet" and lambda.
#'
#'
#'@examples
#'n=20
#'VARknown <- rWishart(1,df=4,Sigma=matrix(c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3))
#'Y0 <- mvrnorm(n=n,mu=rep(0.5,3),Sigma=VARknown[,,1])
#'fitlasso <- lassoglmnet(Y0=Y0,alpha=0.5)
#'
#'@export


lassoglmnet <- function(Y0,nlambda=10,alpha)
{
    p <- dim(Y0)[2]
    n <- dim(Y0)[1]
    coef_fit_glmnet <- array(rep(0,nlambda*p*p),dim=c(nlambda,p,p));
    resid_fit_glmnet <- array(dim=c(nlambda,p,n))
    
    for (i in seq_len(p))
    {
        Y_X=Y0[,-i];
        Y=Y0[,i]
        #glmnet
        glmnetfit <- glmnet(Y_X,Y,intercept=TRUE,nlambda=5,alpha=alpha)
        maxlambda <- glmnetfit$lambda[1]
        minlambda <- glmnetfit$lambda[5]
        
        fit_glmnet <- glmnet(Y_X,Y,intercept=TRUE,
        lambda=seq(from=0,to=maxlambda,length.out=nlambda),alpha=alpha,
        lambda.min.ratio=0.1)
        B_glmnet <- as.matrix(coef(fit_glmnet))
        coef_fit_glmnet[,2:p,i] <- t(B_glmnet)[,2:p] #exclude intercept term
        response <- predict.glmnet(newx=Y_X,fit_glmnet,
        type="response",exact=TRUE)
        #row=obs and col=slambda
        resid_fit_glmnet[,i,] <- t(response-Y);
        lambda <- fit_glmnet$lambda
    }
        return(list(coef_glmnet=coef_fit_glmnet,resid_glmnet=resid_fit_glmnet,
                    lambda=lambda))
}
