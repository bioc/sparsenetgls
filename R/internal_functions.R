## The internal functions for sparsenetgls
## @description Internal sparsenetgls functions
## @details It is not to be called by users.  @type
## others

density <- function(A, p) {
    cn = 0
    A_vals = as.matrix(A)[lower.tri(A)]
    cn = sum(A_vals != 0)
    return(cn)
}


rowVar <- function(x) {
    varx <- apply(x, 1, function(x) var(x, na.rm = TRUE))
    return(varx)
}

colMedian <- function(x) {
    ncol <- dim(x)[2]
    medianx <- apply(x, 2, function(x) median(x, na.rm = TRUE))
    return(medianx)
}


# Distance
poweradj <- function(adj, power) {
    A = adj
    for (i in seq_len(power - 1)) A = A %*% adj
    return(A)
}

convert_to_adj <- function(A, p) {
    A_vals = unlist(lapply(A[lower.tri(A)], function(x) if (x != 
        0) 
        x = 1 else x = 0))  #A is symmetrical
    Am = matrix(nrow = p, ncol = p)
    # set the diagonal term to 0
    ij = cbind(seq_len(p), seq_len(p))
    Am[ij] = 0
    # set the lower and upper triangle
    Am[lower.tri(Am)] = A_vals
    ll = t(Am)
    Am[upper.tri(Am)] = ll[upper.tri(ll)]
    return(Am)
}

add_connect <- function(adjlast, adjnew, p = p) {
    adjnewc <- convert_to_adj(adjnew, p = p)
    diff <- adjnewc - adjlast
    nodeset = c(0, 0)
    diff_vals = diff[lower.tri(diff)]
    adjlast_vals = adjlast[lower.tri(adjlast)]
    
    # index the lower triangle
    p_len <- seq_len(p)
    j <- rep.int(p_len, rev(p_len))
    i <- unlist(lapply(p_len, function(x) p_len[x:p]))
    
    ij = cbind(i, j)
    ij = ij[i != j, ]
    k <- which((diff_vals > 0) & (adjlast_vals != 1))
    
    adjlast_vals[((diff_vals > 0) & (adjlast_vals != 
        1))] = 1
    
    # make it symmetrical
    adjlast[lower.tri(adjlast)] = adjlast_vals
    ll <- t(adjlast)
    adjlast[upper.tri(adjlast)] = ll[upper.tri(ll)]
    
    nodeset <- rbind(nodeset, ij[k, ])
    
    m <- nrow(nodeset)
    if (m > 1) 
        nodeset_select = nodeset[2:m, ] else nodeset_select = "NA"
    
    return(list(adj = adjlast, nodeset = nodeset_select))
}

convert_cov <- function(adj, varY, p = p) {
    varY = as.matrix(adj * varY)
    return(varY)
}

convert_prec <- function(adjm, nlambda, sample_var, 
    Y, p = p) {
    dv <- diag(sample_var)
    lmatrix <- matrix(nrow = p, ncol = p)
    PREC_seq <- lapply(seq(nlambda), function(x) lmatrix)
    I <- Diagonal(p, rep(1, p))
    
    PREC_seq = lapply(adjm, function(x) {
        cory <- x * cor(Y) + I
        return((sqrt(dv)^(-1) * I) %*% ginv(as.matrix(cory)) %*% 
            (sqrt(dv)^(-1) * I))
    })
    
    return(PREC_seq)
}

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


dist_tune <- function(covy, covstart, ndist, p) {
    # Distance tuning
    Aest <- convert_to_adj(covstart, p = p)  #Use the given starting cov matrix
    cov_adj <- array(dim = c(p, p, ndist))
    power = 1
    
    adjnew = Aest
    q <- density(Aest, p = p) + 1  #Giving a safe starting value for q
    
    while (q > 0 & power <= ndist) {
        adjnew <- poweradj(adj = Aest, power = power)
        if (power == 1) {
            cum_adj = Aest
            cov_adj[, , power] <- convert_cov(adj = Aest, 
                varY = covy, p = p)
            +as.matrix(Diagonal(n = p, x = diag(covy)))
        } else {
            cum_connection <- add_connect(adjlast = cum_adj, 
                adjnew = adjnew, p = p)
            # Output cum_connection adj and nodeset
            nodeset <- cum_connection$nodeset
            cum_adj <- cum_connection$adj
            sig_diag <- as.matrix(Diagonal(n = p, x = diag(covy)))
            cov_adj[, , power] <- convert_cov(adj = cum_adj, 
                varY = covy, p = p) + sig_diag
        }
        # new-pairs of nodes added by add_connect
        if (power > 1) {
            q <- nrow(nodeset)
            if (is.null(q)) {
                if (length(nodeset) > 1) {
                q = 1
                nodeset <- matrix(nodeset, nrow = 1, 
                ncol = 2)
                } else q = 0
            }
        }
        power = power + 1
    }
    return(list(cov_adj = cov_adj, power = (power - 
        1)))
}

beta_to_omega <- function(Beta, resid, pathnumber) {
    p <- dim(Beta)[2]
    n <- dim(resid)[1]
    OM_glmnet <- array(rep(0, pathnumber * p * p), 
        dim = c(pathnumber, p, p))
    OMEGA_glmnet <- array(dim = c(p, p, pathnumber))
    
    resid_vals <- matrix(as.vector(resid), nrow = n, 
        ncol = pathnumber * p)
    
    # resid is an array with c(n,p,nlambda)
    dm <- matrix(apply(resid_vals, 2, function(x) var(x, 
        na.rm = TRUE)), ncol = pathnumber, nrow = p)
    
    # Array operation beta
    # array->coef_fit_glmnet[beta_2:beta_p,nodes,lambda]
    
    for (l in seq_len(pathnumber)) {
        Beta_vals_L = Beta[, , l][lower.tri(Beta[, 
            , l])]
        Beta_vals_U = t(Beta[, , l])[lower.tri(t(Beta[, 
            , l]))]
        j = seq_len(p - 1)
        k = seq(p - 1, 1)  #rep times 
        dm_U = unlist(lapply(j, function(x) rep(dm[seq_len(p - 
            1), l][x], each = k[x])))
        j = seq_len(p - 1) + 1
        dm_L = unlist(lapply(j, function(x) rep(dm[seq(x, 
            p), l], 1)))
        
        # assign the lower triangle of OMEGA
        sign_same = (sign(Beta_vals_L) == sign(Beta_vals_U))
        OM_glmnet_vals = vector(mode = "numeric", length = length(Beta_vals_L))
        OM_glmnet_vals[sign_same] = (-1) * sign(Beta_vals_L[sign_same]) * 
            sqrt((Beta_vals_L[sign_same]/dm_L[sign_same]) * 
                (Beta_vals_U[sign_same]/dm_U[sign_same]))
        OM_glmnet_vals[!sign_same] = 0
        
        OMEGA_glmnet[, , l][lower.tri(OMEGA_glmnet[, 
            , l])] = OM_glmnet_vals
        ll = t(OMEGA_glmnet[, , l])
        OMEGA_glmnet[, , l][upper.tri(OMEGA_glmnet[, 
            , l])] = ll[upper.tri(ll)]
        diag(OMEGA_glmnet[, , l]) = 1/dm[, l]
    }
    
    return(OMEGAMATRIX = OMEGA_glmnet)
}


prec_select <- function(prec, precstart, p = p) {
    # selecting precision matrix terms
    
    Aest <- convert_to_adj(precstart, p = p)
    # Use the given starting prec matrix
    prec_selected <- prec * Aest + as.matrix(Diagonal(n = p, 
        diag(prec)))
    
    return(prec_selected)
}

