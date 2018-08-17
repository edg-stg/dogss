iter_pars <- function(pars) {
    
    X <- pars$X
    Y <- pars$Y
    G <- pars$G
    sigma2_0 <- pars$sigma2_0
    sigma2_slab <- pars$sigma2_slab
    M <- pars$M
    N <- pars$N
    nbyG <- pars$nbyG
    groups <- pars$groups
    damping <- pars$damping
    
    if (M != length(groups)) {
        # if length(groups)==M, we don't have any groups, so we don't update
        # anything for f3tilde and ignore f4tilde update for f3tilde first... find
        # pars of marginal q3slash distribution
        rho3_slash <- pars$rho[G] - pars$rho3
        r3_slash <- pars$r - pars$r3  # = pars$r2 !!!
        rho3_new <- log(0.5) + logsum1pexp(r3_slash)  # log(0.5) + log(1/(1-p3_slash))  # -log(1-p3_slash) + log(p3_slash*p0 + (1-p3_slash)*(1-p0))
        r3_new <- -logsum1pexp(-rho3_slash + log(2))  # log(pi3_slash*p0/(1-pi3_slash*p0)) = logit(pi3_slash*p0)
        
        # damping!
        pars$rho3 <- damping * rho3_new + (1 - damping) * pars$rho3
        pars$r3 <- damping * r3_new + (1 - damping) * pars$r3  # eigtl p3^alpha*p3old^(1-alpha)
        
        # update parameter p and pi (resp. r and rho) of posterior:
        r_new <- pars$r2 + pars$r3
        rho_new <- pars$rho4 + sapply(groups, function(g) sum(pars$rho3[G == 
            g]))
        pars$r <- r_new
        pars$rho <- rho_new
    }
    
    # find pars of the marginal q-slash-2-distributions
    V2_slash_new <- 1/(1/diag(pars$V) - pars$V2_inv)
    V2_slash_new[!is.finite(V2_slash_new)] <- 1e+06  # sometimes variance too large
    ## if there are terms with negative variance V2_slash, we do not perform the
    ## update for these terms (resp. doing it with the old V2_slash?):
    indices <- which(V2_slash_new > 0)
    V2_slash <- pars$V2  #_slash
    V2_slash[indices] <- V2_slash_new[indices]
    
    m2_slash <- V2_slash * (1/diag(pars$V) * pars$m - pars$V2_inv_m2)
    
    r2_slash <- pars$r - pars$r2
    
    # derive pars of the q-new-distribution: not needed (but p_aux)! derive
    # updated pars of f2tilde directly:
    
    r2_new <- 0.5 * (log(V2_slash) - log(V2_slash + sigma2_slab) + m2_slash^2 * 
        (1/V2_slash - 1/(V2_slash + sigma2_slab)))
    p_aux <- sigmoid(r2_new + r2_slash)  # , pars$sigmoid_function
    
    a <- p_aux * m2_slash/(V2_slash + sigma2_slab) + (1 - p_aux) * m2_slash/V2_slash
    b <- p_aux * (m2_slash^2 - V2_slash - sigma2_slab)/(V2_slash + sigma2_slab)^2 + 
        (1 - p_aux) * (m2_slash^2 - V2_slash)/V2_slash^2
    V2_new <- 1/(a^2 - b) - V2_slash
    m2_new <- m2_slash - a * (V2_new + V2_slash)
    
    # get rid of variances V2 equal to zero or negative:
    V2_new[V2_new == 0] <- 1e-10
    V2_new[V2_new < 0] <- 100  # why 100? see Hernandezlobato 2010
    
    # damping!
    pars$V2_inv <- damping * 1/V2_new + (1 - damping) * pars$V2_inv
    pars$V2 <- 1/pars$V2_inv
    pars$V2_inv_m2 <- damping * 1/V2_new * m2_new + (1 - damping) * pars$V2_inv_m2
    pars$r2 <- damping * r2_new + (1 - damping) * pars$r2  # eigtl p2^alpha*p2old^(1-alpha)
    
    pars$damping <- damping * 0.99
    
    # update pars m and V and p of posterior:
    nu <- pars$V1_inv_m1 + pars$V2_inv_m2
    
    if (M > N) {
        # more features than observations: Woodbury trick Lambda <-
        # diag(as.vector(pars$V2), nrow=M) Sigma2D <- diag(sigma2_0, N)
        pars$V <- CppEigenWoodbury(pars$V2, X, rep(sigma2_0, N))  # upper triangular matrix, thats enough
        pars$m <- CppEigenVtimesNu(pars$V, nu)  # take advantage of upper triangular/symmetric matrix
        # pars$V[lower.tri(pars$V)] <- t(pars$V)[lower.tri(pars$V)] rm(Lambda,
        # Sigma2D) XLambda <- t(t(X)*as.vector(pars$V2)) # that's the same like X
        # %*% Lambda, but much more efficient sigma2pXLambdaXt <- diag(sigma2_0, N)
        # + CppEigenMultiplication(XLambda, t(X)) # XLambda %*% t(X) #
        # tcrossprod(XLambda, X) # Inv <- CppEigenInv(sigma2pXLambdaXt) #
        # solve(sigma2pXLambdaXt, tol = 1e-100) rm(sigma2pXLambdaXt) LInv <- Inv
        # LInv[upper.tri(LInv)] <- 0 diag(LInv) <- diag(LInv)/2 # t(XLambda) %*%
        # LInv %*% XLambda pars$V <- Lambda - CppEigenMultiplicationThree(XLambda,
        # LInv) # CppMultiplicationThree(t(XLambda), Inv, XLambda) #
        # tcrossprod(crossprod(XLambda, Inv), t(XLambda)) # rm(Lambda, XLambda,
        # Inv) C <- tryCatch({chol(sigma2pXLambdaXt)}, error = function(e)
        # {return(NULL)}) Inv <- 'if'(is.null(C), solve(sigma2pXLambdaXt),
        # chol2inv(C)) #for solve: , tol = 1e-100 if (is.null(C)) cat('cholesky
        # didn't work!\n') V <- Lambda - t(XLambda) %*% Inv %*% XLambda more
        # observations N than features M: calculate V directly
    } else {
        pars$V <- CppEigenInv(pars$V1_inv + diag(1/as.vector(pars$V2), nrow = M))  # solve(pars$V1_inv+diag(1/as.vector(pars$V2), nrow=M), tol=1e-100)
        pars$m <- as.vector(pars$V %*% nu)
        # C <- tryCatch({chol(V1_inv+diag(1/V2))}, error = function(e)
        # {return(NULL)}) Inv <- 'if'(is.null(C), solve(V1_inv+diag(1/V2)),
        # chol2inv(C)) #for solve: , tol = 1e-100 if (is.null(C)) cat('cholesky
        # didn't work!\n') V <- Inv
    }
    
    pars$r <- pars$r2 + pars$r3
    
    if (any(pars$r > 0)) {
        pars$predicterror <- sum((Y - X[, pars$r > 0, drop = FALSE] %*% pars$m[pars$r > 
            0])^2)
    } else {
        pars$predicterror <- sum(Y^2)
    }
    
    
    pars$iterations <- pars$iterations + 1
    
    return(pars)
    
}
