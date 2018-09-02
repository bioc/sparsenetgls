#'The sparsenetgls() function
#'@description The sparsenetgls functin is a combination of the graph structure
#'learning and generalized least square regression.
#'It is designed for multivariate regression uses penalized and/or regularised
#'approach to deriving the precision and covariance Matrix of the multivariate
#'Gaussian distributed responses. Gaussian Graphical model is used to learn the
#'structure of the graph and construct the precision and covariance matrix.
#'Generalized least squared regression is used to derive the sandwich
#'estimation of variances-covariance matrix for the regression coefficients
#'of the explanatory variables, conditional on the solutions of the precision
#'and covariance matrix.
#'
#'@import Matrix MASS huge parcor
#'@importFrom methods new
#'@importFrom stats coef cor cov median sd var
#'@importFrom utils stack
#'@importFrom graphics legend lines par plot points segments
#'
#'@param responsedata It is a data matrix of multivariate normal distributed
#'response variables. Each row represents one observation sample and each
#'column represents one variable.
#'
#'@param predictdata It is a data matrix of explanatory variables and has the
#'same number of rows as the response data.
#'
#'@param nlambda It is an interger recording the number of lambda value used
#'in the penalized path for estimating the precision matrix.
#'The default value is 10.
#'
#'@param ndist It is an interger recording the number of distant value used in
#'the penalized path for estimating the covariance matrix.
#'The default value is 5.
#'
#'@param method It is the option parameter for selecting the penalized method
#'to derive the precision matrix in the calculation of the sandwich estimator
#'of regression coefficients and their variance-covariance
#'matrix. The options are 'glasso', 'lasso','elastic', and 'mb'. 'glasso' use
#'the graphical lasso method documented in Yuan and lin (2007) and Friedman,
#'Hastie et al (2007). It used the imported function from
#'R package 'huge'. 'lasso' use the penalized liner regression among the
#'response variables (Y[,j]~Y[,1]+...Y[,j-1],Y[,j+1] +...Y[,p]) to estimate
#'the precision matrix.  'elastic' uses the enet-regularized linear
#'regression among the response variables to estimate the precision matrix.
#'Both of these methods utilize the coordinate descending algorithm documentd
#'in Friedman, J., Hastie, T. and Tibshirani, R. (2008) and use the
#'imported function from R package 'glmnet'.  'mb' use the Meinshausen and
#'Buhlmann penalized linear regression and the neighbourhood selection with
#'the lasso approach (2006) to select the covariance terms and derive the
#'corresponding precision matrix ; It uses the imported function from 'huge'
#'in function sparsenetgls().
#'
#'@param lambda.min.ratio It is the default parameter set in function huge()
#'in the package 'huge'. Quoted from huge(),
#'it is the minial value of lambda, being a fraction of the uppperbound (MAX)
#'of the regularization/thresholding parameter
#'that makes all the estimates equal to 0. The default value is 0.001.
#'It is only applicable when 'glasso' and 'mb' method is used.
#'
#'@return  Return the list of regression results including the regression
#'coefficients, array of variance-covariance matrix
#'for different lambda and distance values, lambda and distance (power) values,
#'bic and aic for model fitting, and the
#'list of precision matrices on the penalized path.
#'
#'@examples
#'ndox=5; p=3; n=1000
#'VARknown <- rWishart(1, df=4, Sigma=matrix(c(1,0,0,0,1,0,0,0,1),
#'nrow=3,ncol=3))
#'normc <- mvrnorm(n=n,mu=rep(0,p),Sigma=VARknown[,,1])
#'Y0=normc
#'##u-beta
#'u <- rep(1,ndox)
#'X <- mvrnorm(n=n,mu=rep(0,ndox),Sigma=Diagonal(ndox,rep(1,ndox)))
#'X00 <- scale(X,center=TRUE,scale=TRUE)
#'X0 <- cbind(rep(1,n),X00)
#'#Add predictors of simulated CNA
#'abundance1 <- scale(Y0,center=TRUE,scale=TRUE)+as.vector(X00%*%as.matrix(u))
#'
#'##sparsenetgls()
#'fitgls <- sparsenetgls(responsedata=abundance1,predictdata=X00,
#'nlambda=5,ndist=2,method='glasso')
#'nlambda=5
#'##rescale regression coefficients from sparsenetgls
#'#betagls <- matrix(nrow=nlambda, ncol=ndox+1)
#'#for (i in seq_len(nlambda))   
#'#betagls[i,] <- convertbeta(Y=abundance1, X=X00, q=ndox+1,
#'#beta0=fitgls$beta[,i])$betaconv
#'@export
sparsenetgls <- function(responsedata, predictdata, 
    nlambda = 10, ndist = 5, method = c("lasso", "glasso", 
        "elastic", "mb"), lambda.min.ratio = 1e-05) {
    
    method <- match.arg(method)
    colnameprotein <- colnames(responsedata)
    # Add predictors(subset Y)
    p <- dim(responsedata)[2]
    ndox <- dim(predictdata)[2]
    n <- dim(responsedata)[1]
    Y1 <- scale(responsedata)
    X0 <- cbind(rep(1, n), scale(predictdata))
    
    # Stack Y and X
    nodes <- seq_len(p)
    units <- seq_len(n)
    colnames(Y1) <- nodes
    response <- data.frame(Y1)
    responsedata <- stack(response)
    index <- rep(seq_len(p), n)
    responsedata$ind <- sort(index)
    stunits <- rep(units, p)
    X <- apply(X0, 2, rep, p)
    colnames(X) <- colnames(X0)
    Y <- responsedata$values
    
    beta0 <- solve(t(X) %*% X) %*% crossprod(X, Y)
    covy <- cov(Y1 - as.vector(crossprod(t(X0), beta0)))
    
    D <- as.matrix(Diagonal(x = diag(covy)))
    
    beta <- array(dim = c(ndox + 1, nlambda))
    covBeta <- array(dim = c(ndox + 1, ndox + 1, ndist, 
        nlambda))
    
    # This serie uses the sample data PREC and
    # COV(Y1-X%*%BETA)
    if (method == "glasso") {
        penalized_seq <- huge(covy, lambda.min.ratio = lambda.min.ratio, 
            nlambda = nlambda, scr = FALSE, cov.output = TRUE, 
            method = "glasso")
        PREC_seq = penalized_seq$icov
        lambdav = penalized_seq$lambda
        
    } else if (method == "mb") {
        penalized_seq <- huge(Y1, lambda.min.ratio = lambda.min.ratio, 
            nlambda = nlambda, scr = TRUE, sym = "and", 
            method = "mb")
        PREC_seq <- convert_prec(adjm = penalized_seq$path, 
            nlambda = nlambda, sample_var = covy, Y = Y1, 
            p = p)
        lambdav = penalized_seq$lambda
        
    } else if (method == "lasso") {
        fitlasso <- lassoglmnet(Y1, nlambda = nlambda, 
            alpha = 1)
        lambdav = fitlasso$lambda
        
        PREC_ar <- beta_to_omega(Beta = fitlasso$coef_glmnet, 
            resid = fitlasso$resid_glmnet, pathnumber = nlambda)
        PREC_seq <- lapply(seq(nlambda), function(x) PREC_ar[, 
            , x])
        rm(PREC_ar)
        
    } else if (method == "elastic") {
        fitlasso_enet <- lassoglmnet(Y1, nlambda = nlambda, 
            alpha = 0.5)
        lambdav = fitlasso_enet$lambda
        PREC_ar <- beta_to_omega(Beta = fitlasso_enet$coef_glmnet, 
            resid = fitlasso_enet$resid_glmnet, pathnumber = nlambda)
        PREC_seq <- lapply(seq(nlambda), function(x) PREC_ar[, 
            , x])
        rm(PREC_ar)
    }
    
    # Derive the path sandwitch parameter and its
    # variance using the chosen cov of tuning parameter
    # d
    dist_adj <- dist_tune(covy = covy, covstart = PREC_seq[[round(nlambda/2, 
        0)]], ndist = ndist, p = p)
    cov_adj <- dist_adj$cov_adj
    
    B = array(dim = c(ndox + 1, ndox + 1, ndist, nlambda))
    aic = bic = array(dim = c(dist_adj$power, nlambda))
    
    for (lambda in seq_len(nlambda)) {
        power = 1
        while (power <= dist_adj$power) {
        cov_adj_kron <- as.matrix(bdiag(rep(list(cov_adj[,, power]), n)))
        V_1_kron <- as.matrix(bdiag(rep(list(as.matrix(PREC_seq[[lambda]])),
        n)))
        B[, , power, lambda] <- t(X) %*% V_1_kron %*% cov_adj_kron %*%
        (V_1_kron %*% X)
        power = power + 1
        }
    }
    
    for (lambda in seq_len(nlambda)) {
        for (s in seq_len(dist_adj$power)) {
        V_1_kron <- as.matrix(bdiag(rep(list(as.matrix(PREC_seq[[lambda]])),n)))
        L <- solve(t(X) %*% V_1_kron %*% X)
        H <- t(X) %*% as.matrix(V_1_kron)
        beta[, lambda] <- L %*% (H %*% Y)
        covBeta[, , s, lambda] <- L %*% B[, , s,lambda] %*% L
        bic[s, lambda] <- -Matrix::determinant(PREC_seq[[lambda]])$modulus[1]
        +t(Y - X %*% beta[, lambda]) %*% V_1_kron %*%(Y - X %*% beta[, lambda])
            +density(PREC_seq[[lambda]], p) * log(n)/n
        aic[s, lambda] <- -Matrix::determinant(PREC_seq[[lambda]])$modulus[1]
            +t(Y - X %*% beta[, lambda]) %*% V_1_kron %*% 
            (Y - X %*% beta[, lambda])
            +density(PREC_seq[[lambda]], p)/n
        }
    }
    return(list(beta = beta, covBeta = covBeta[, , 
        seq_len(dist_adj$power), ], lambda = lambdav, 
        power = power, bic = bic, aic = aic, PREC_seq = PREC_seq))
}










