#'The internal functions for sparsenetgls
#'@description  Internal sparsenetgls functions
#'@details It is not to be called by users.
#'@param A A matrix for counting the edges in density() function
#'@param p Number of dimension in density() function
#'@return Return half no of the links of Matrix A
#'@import methods Matrix MASS
#'@importFrom grDevices dev.new dev.off
#'@importFrom graphics legend lines par plot points segments
#'@importFrom methods new
#'@importFrom stats coef cor cov median sd var
#'@importFrom utils stack


density <- function(A,p)
{cn=0
for (i in seq_len(p))
    for (j in seq_len(i))
    {if ((i!=j) & (A[i,j]!=0)) cn=cn+1}
return(cn)
}

rowVar <- function(x)
{
    nrow <- dim(x)[1]
    varx <- vector("numeric",length=nrow)
    for (j in seq_len(nrow)) varx[j]=var(x[j,],na.rm=TRUE)
return(varx)
}

colMedian <- function(x)
{
    ncol <- dim(x)[2]
    medianx <- vector("numeric",length=ncol)
    for (j in seq_len(ncol)) medianx[j] <- median(x[,j],na.rm=TRUE)
return(medianx)
}

#Distance
poweradj <- function(adj,power){
A=adj;
for ( i in seq_len(power-1)) A=A%*%adj
return(A);
}

convert_to_adj <- function(A,p){
for (i in seq_len(p))
    for (j in seq_len(i))    #Including the diagnal term if from 1 to i
    {if ((i!=j) & (A[i,j]!=0)) A[i,j]=A[j,i]=1 else if (i==j) A[i,j]=0}
return(A)
}

add_connect <- function(adjlast,adjnew,p=p){
    adjnewc <- convert_to_adj(adjnew,p=p)
    diff=adjnewc-adjlast
    nodeset=c(0,0)
    for ( i in seq_len(p))
        for ( j in seq_len(i))
            if ((diff[i,j]>0) & (adjlast[i,j]!=1))
        {adjlast[i,j]=adjlast[j,i]=1
        nodeset=rbind(nodeset,c(i,j))}    #lower triangle
    m <- nrow(nodeset)
if (!is.null(m)) nodeset_select=nodeset[2:m,] else nodeset_select='NA'
    return(list(adj=adjlast,nodeset=nodeset_select))
}

convert_cov <- function(adj,varY,p=p)
{
for ( i in seq_len(p))
    for ( j in seq_len(i))
        if (adj[i,j]==0) varY[i,j]=varY[j,i]=0
        return(varY)
}


convertbeta <- function(X,Y,q,beta0) {
    betaconv=beta0
    meanX <- colMeans(X,na.rm=TRUE)
    meanY <- mean(Y,na.rm=TRUE)
    sdX <- sqrt(rowVar(t(X)))
    sdY <- apply(Y,2,sd,na.rm=TRUE)
    conv_ratio <- (mean(sdY)*sdX^(-1))
    betaconv_int<- -conv_ratio*meanX*beta0[2:q]+beta0[1]*mean(sdY)+meanY
    betaconv[1] <- mean(betaconv_int,na.rm=TRUE)
    betaconv[2:q] <- conv_ratio*beta0[2:q]
    return(list(betaconv=betaconv,betaconv_int=betaconv_int))}
    
convert_prec <- function(adjm,nlambda,sample_var,Y,p=p) {
dv <- diag(sample_var)
lmatrix <- matrix(nrow=p,ncol=p)
PREC_seq <- lapply(seq(nlambda),function(x) lmatrix)
I <- Diagonal(p,rep(1,p))

for (i in seq_len(nlambda)) {
cory <- adjm[[i]]*cor(Y)+I
PREC_seq[[i]] <- (sqrt(dv)^(-1)*I)%*%ginv(as.matrix(cory))%*%(sqrt(dv)^(-1)*I)
    }
return (PREC_seq)
}


dist_tune <- function(covy,covstart,ndist,p)
{
#Distance tuning
    Aest <- convert_to_adj(covstart,p=p)   #Use the given starting cov matrix
    cov_adj <- array(dim=c(p,p,ndist))
    power=1
    
    adjnew=Aest
    q <- density(Aest,p=p)+1               #Giving a safe starting value for q
    
    while (q>0 & power<=ndist)
    {
    adjnew <- poweradj(adj=Aest,power=power)
    if (power==1) {
        cum_adj=Aest
        cov_adj[,,power] <- convert_cov(adj=Aest,varY=covy,p=p)
        +as.matrix(Diagonal(n=p,x=diag(covy)))}  else
    {
    cum_connection <- add_connect(adjlast=cum_adj,adjnew=adjnew,p=p)
                                #Output cum_connection adj and nodeset
        nodeset <- cum_connection$nodeset
        cum_adj <- cum_connection$adj
        sig_diag <- as.matrix(Diagonal(n=p,x=diag(covy)))
        cov_adj[,,power] <- convert_cov(adj=cum_adj,varY=covy,p=p)+sig_diag
    }
    #new-pairs of nodes added by add_connect
        if (power>1)
        {q <- nrow(nodeset)
        if (is.null(q)) {if (length(nodeset)>1) {
            q=1;nodeset <- matrix(nodeset,nrow=1,ncol=2)} else q=0}
        }
        power=power+1;
    }
return(list(cov_adj=cov_adj,power=(power-1)))
}


beta_to_omega <- function(Beta,resid,pathnumber)
{    
    p <- dim(Beta)[2]
    OM_glmnet <- array(rep(0,pathnumber*p*p),dim=c(pathnumber,p,p))
    OMEGA_glmnet <- array(dim=c(p,p,pathnumber))
    dm <- matrix(nrow=pathnumber,ncol=p)
    
    for (i in seq_len(p))
        for ( l in seq_len(pathnumber)) {dm[l,i] <- var(resid[l,i,])}

    for ( l in seq_len(pathnumber))
        for (j in seq_len(p-1))
            #j represents node and column stores coef of j node
            for ( i in 2:p)           #only uses left lower triangle
            if (i>j) {
                if (sign(Beta[l,i,j])==sign(Beta[l,(j+1),i])) {
                pcor=sqrt((Beta[l,i,j]/dm[l,j])*(Beta[l,(j+1),i]/dm[l,i]))
                OM_glmnet[l,i,j]=(-1)*sign(Beta[l,i,j])*pcor
                } else 
                OM_glmnet[l,i,j]=0
                }
            #sparsed solution will not give the symmetrical OMEGA

    for ( l in seq_len(pathnumber))
            for (j in seq_len(p))
                for (i in seq_len(p))
            {if (i>j)  
            {OMEGA_glmnet[i,j,l]=OMEGA_glmnet[j,i,l]=OM_glmnet[l,i,j]}
                if (i==j) {OMEGA_glmnet[i,j,l]=1/dm[l,i]}
                }
    return(OMEGAMATRIX=OMEGA_glmnet)
    }


prec_select <- function(prec,precstart,p=p)
{
#selecting precision matrix terms
    
    Aest <- convert_to_adj(precstart,p=p)
#Use the given starting prec matrix
    prec_selected <- prec*Aest+as.matrix(Diagonal(n=p,diag(prec)))
    
return(prec_selected)
}

#'Fast version of Matrix :: .bdiag() -- for the case of *many* (k x k)
#'matrices:
#'@import methods Matrix MASS
#'@param lmat list(<mat1>, <mat2>, ....., <mat_N>) where each mat_j is
#'a k x k 'matrix'
#'@return a sparse (N*k x N*k) matrix of class
#'\code{"\linkS4class{dgCMatrix}"}.
#'

bdiag_m <- function(lmat) {
## Copyright (C) 2016 Martin Maechler, ETH Zurich
    
if(!length(lmat)) return(new("dgCMatrix"))
    stopifnot(is.list(lmat), is.matrix(lmat[[1]]),
            (k <- (d <- dim(lmat[[1]]))[1]) == d[2], # k x k
            all(vapply(lmat, dim, integer(2)) == k)) # all of them
    N <- length(lmat)
    if(N * k > .Machine$integer.max)
    stop("resulting matrix too large; would be M x M, with M=", N*k)
    M <- as.integer(N * k)
## result: an M x M matrix
    new("dgCMatrix", Dim = c(M,M),
    ## 'i :' maybe there's a faster way (w/o matrix indexing), but elegant?
    i <- as.vector(matrix(0L:(M-1L), nrow=k)[, rep(seq_len(N), each=k)]),
    p <- k * 0L:M,
    x <- as.double(unlist(lmat, recursive=FALSE, use.names=FALSE)))
}
