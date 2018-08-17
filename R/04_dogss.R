#' @export dogss

dogss <- function(X, Y, G = NULL, sigma_0 = 1, sigma_slab = 2, damping = 0.9,
    standardize = FALSE, intercept = FALSE, timeseries = FALSE, maxdelay = 1,
    tol = 1e-05, iter.max = 100) {

    parameters <- init_pars(X, Y, G, sigma_0, sigma_slab, damping, standardize,
        intercept, timeseries, maxdelay)

    convergence <- FALSE
    while (!convergence) {
        m_old <- parameters$m
        predicterror_old <- parameters$predicterror
        parameters <- iter_pars(parameters)
        convergence <- converged(m_old, parameters$m, predicterror_old, parameters$predicterror,
            parameters$iterations, tol, iter.max)
    }

    if (parameters$M != length(parameters$groups)) {
        parameters$p_final <- sigmoid(sapply(parameters$G, function(g) parameters$rho[parameters$groups ==
            g])) * sigmoid(parameters$r)  # , sigmoid_function
    } else {
        # if length(groups)==M, we don't have any groups, corresponds to simple
        # spikeandslab
        parameters$p_final <- sigmoid(parameters$r)
    }
    names(parameters$p_final) <- names(parameters$r)

    return(list(p_final = parameters$p_final, m = parameters$m, p_features = sigmoid(parameters$r),
        p_groups = sigmoid(parameters$rho), sigma2_slab = parameters$sigma2_slab,
        sigma2_0 = parameters$sigma2_0, X = parameters$X, Y = parameters$Y,
        G = parameters$G, groups = parameters$groups, intercept = parameters$intercept,
        timeseries = parameters$timeseries, maxdelay = parameters$maxdelay,
        iterations = parameters$iterations, prederror = parameters$predicterror))

}
