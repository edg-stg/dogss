init_pars <- function(X, Y, G = NULL, sigma_0 = 1, sigma_slab = 2, damping = 0.9, 
    standardize, intercept, timeseries, maxdelay) {
    # , sigmoid_function='exact'
    
    N <- length(Y)
    
    if (is.null(dim(X))) 
        X <- as.matrix(X, ncol = 1)
    
    if (standardize) 
        X <- scale(X)
    
    if (timeseries) {
        nvar <- ncol(X)
        X <- expandMatrix(X, maxdelay)
        Y <- Y[c((maxdelay + 1):N)]
        N <- N - maxdelay
        G <- rep(1:nvar, each = maxdelay)
    }
    
    if (intercept) {
        X <- cbind(rep(1, dim(X)[1]), X)
        if (is.null(G)) {
            G <- c("intercept", 1:(dim(X)[2] - 1))
        } else {
            G <- c("intercept", G)
        }
    }
    
    if (is.null(G)) 
        G <- 1:dim(X)[2]
    
    M <- length(G)  # number of features (=columns of X or length of G)
    sigma2_0 <- sigma_0 * sigma_0
    sigma2_slab <- sigma_slab * sigma_slab
    
    G <- as.character(G)
    
    groups <- unique(G)
    
    nbyG <- sapply(groups, function(x) sum(G == x))  # number of variables in every group
    names(nbyG) <- groups
    
    p0 <- rep(0.5, M)  # rep(1/M, M) # no prior: chance of 0.5 to include variable
    names(p0) <- G
    pi0 <- rep(0.5, length(groups))  # rep(1/length(groups), length(groups)) # no prior: chance of 0.5 to include group
    names(pi0) <- groups
    # probabilities transformed for stability:
    r0 <- logit(p0)
    rho0 <- logit(pi0)
    rho <- rho0
    # rho2 <- rho0
    rho3 <- r0  # !!! there's a 'group' probability for every feature in the marginal f3tilde-distribution
    rho4 <- rho0
    
    # these parameters won't change by EP:
    V1_inv <- 1/sigma2_0 * t(X) %*% X
    V1_inv_m1 <- as.vector(1/sigma2_0 * t(X) %*% Y)
    
    # initialize the remaining parameters:
    V2 <- p0[1:M] * sigma2_slab  # V2 and V2_inv are diagonal matrices, stored as vectors # why p_0*sigma_slab?
    V2_inv <- 1/V2
    V2_inv_m2 <- rep(0, M)
    r2 <- r0
    r3 <- r0
    # names(r3) <- G
    Lambda <- diag(V2, nrow = M)
    nu <- V1_inv_m1 + V2_inv_m2
    
    if (M > N) {
        # more features than observations: Woodbury trick
        XLambda <- t(t(X) * V2)  # that's the same like X %*% Lambda, but much more efficient
        sigma2pXLambdaXt <- diag(sigma2_0, N) + XLambda %*% t(X)
        Inv <- solve(sigma2pXLambdaXt, tol = 1e-100)
        V <- Lambda - t(XLambda) %*% Inv %*% XLambda
        # C <- tryCatch({chol(sigma2pXLambdaXt)}, error = function(e)
        # {return(NULL)}) Inv <- 'if'(is.null(C), solve(sigma2pXLambdaXt),
        # chol2inv(C)) #for solve: , tol = 1e-100 if (is.null(C)) cat('cholesky
        # didn't work!\n') V <- Lambda - t(XLambda) %*% Inv %*% XLambda more
        # observations N than features M: calculate V directly
    } else {
        Inv <- solve(V1_inv + diag(1/V2, nrow = M), tol = 1e-100)
        V <- Inv
        # C <- tryCatch({chol(V1_inv+diag(1/V2))}, error = function(e)
        # {return(NULL)}) Inv <- 'if'(is.null(C), solve(V1_inv+diag(1/V2)),
        # chol2inv(C)) #for solve: , tol = 1e-100 if (is.null(C)) cat('cholesky
        # didn't work!\n') V <- Inv
    }
    
    m <- as.vector(V %*% nu)
    
    r <- r0
    
    predicterror <- sum(Y^2)
    
    # V2_slash <- V2 V2_slash_new <- 1/(1/diag(V) - 1/V2) #### schon im ersten
    # Schritt negativ???  indices <- which(V2_slash_new > 0) V2_slash[indices]
    # <- V2_slash_new[indices]
    
    pars <- list(V1_inv = V1_inv, V1_inv_m1 = V1_inv_m1, p0 = p0, rho4 = rho4, 
        V2 = V2, V2_inv = V2_inv, V2_inv_m2 = V2_inv_m2, r2 = r2, r3 = r3, 
        rho3 = rho3, V = V, m = m, r = r, rho = rho, predicterror = predicterror, 
        sigma2_slab = sigma2_slab, sigma2_0 = sigma2_0, X = X, Y = Y, M = M, 
        N = N, G = G, groups = groups, damping = damping, intercept = intercept, 
        timeseries = timeseries, maxdelay = maxdelay, iterations = 0)  # V2 is needed for Woodbury formula # , sigmoid_function=sigmoid_function
    
    return(pars)
    
}
