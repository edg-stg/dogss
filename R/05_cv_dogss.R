#' @export cv_dogss

cv_dogss <- function(X, Y, G = NULL, nfolds = 10, sigma_0 = 1, sigma_slab = 2,
    damping = 0.9, standardize = FALSE, intercept = FALSE, timeseries = FALSE,
    maxdelay = 1, tol = 1e-05, iter.max = 100) {
    if (timeseries)
        cat("timeseries=TRUE doesn't make sense for cross-validation!\n")

    nobs <- length(Y)
    foldids <- sample(rep(seq(nfolds), length = nobs))
    outlist <- as.list(seq(nfolds))
    probseq <- seq(0.9, 0, -0.1)
    predmat <- matrix(NA, length(Y), length(probseq))
    ms <- matrix(NA, dim(X)[2], nfolds)
    ps <- ms

    for (i in seq(nfolds)) {
        which <- foldids == i
        Y_sub <- Y[!which]
        X_sub <- X[!which, , drop = FALSE]
        outlist[[i]] = dogss(X = X_sub, Y = Y_sub, G = G, sigma_0 = sigma_0,
            sigma_slab = sigma_slab, damping = damping, standardize = standardize,
            intercept = intercept, timeseries = timeseries, maxdelay = maxdelay,
            tol = tol, iter.max = iter.max)
        ms[, i] <- outlist[[i]]$m
        ps[, i] <- outlist[[i]]$p_final
        preds <- matrix(NA, nrow = sum(which), ncol = length(probseq))
        for (p in seq(probseq)) {
            preds[, p] <- X[which, , drop = FALSE] %*% ifelse(ps[, i] >= probseq[p],
                ms[, i], 0)
        }
        predmat[which, ] <- preds
    }

    cvraw <- (Y - predmat)^2
    cvm <- apply(cvraw, 2, mean, na.rm = TRUE)
    cvsd <- sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, mean, na.rm = TRUE)/(nobs -
        1))

    which_cvm_min <- which.min(cvm)
    which_1se <- min(which(cvm <= min(cvm) + cvsd[which_cvm_min]))

    m_final <- apply(ms, 1, mean)
    p_final <- apply(ps, 1, mean)
    m_cv <- ifelse(p_final >= probseq[which_1se], m_final, 0)
    out <- list(m_cv = m_cv, p_final = p_final, m_final = m_final)

    return(out)

}
