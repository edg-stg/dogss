sigmoid <- function(x) {
    # if (type=='exact') {
    y <- double(length(x))
    logic <- (x >= 0)
    # yes <- 1/(1 + exp(-x))
    z <- exp(x[!logic])  # z <- exp(x)
    # no <- z/(1+z)
    y[logic] <- 1/(1 + exp(-x[logic]))  # yes[logic]
    y[!logic] <- z/(1 + z)  # no[!logic]
    return(y)
    # } else if (type=='poly') {# approximation with polynomial kind of y <-
    # double(length(x)) logic <- (x>=16) y[logic] <- 1 tmp <-
    # 1-0.0625*x[!logic] tmp <- tmp*tmp tmp <- tmp*tmp tmp <- tmp*tmp tmp <-
    # tmp*tmp y[!logic] <- 1/(1+tmp) return(y) # } else {# crude approximation
    # with abs return(0.5*(x/(1+abs(x)))+0.5) }
}

# Sigmoid(float signal) (sigmoid(4x)!!!!)  { if (signal >= 4f) return 1f;
# float tmp = 1f - 0.25f * signal; tmp *= tmp; tmp *= tmp; tmp *= tmp; tmp
# *= tmp; return 1f / (1f + tmp); }

logit <- function(x) log(x/(1 - x))

expandMatrix <- function(data, maxdelay) {
    # build a time series matrix
    G <- ncol(data)
    Times <- nrow(data)
    Tdelayed <- Times - maxdelay
    X <- matrix(0, nrow = Tdelayed, ncol = G * maxdelay)
    for (j in 1:maxdelay) {
        for (j2 in 1:G) {
            X[, (j2 - 1) * maxdelay + j] <- data[c((maxdelay - j + 1):(Times - 
                j)), j2]
        }
    }
    return(X)
}

logsumexp <- function(v) {
    # v is a vector if there is at least one element that's not '-Inf', we can
    # do the calculation
    if (sum(is.finite(v)) > 0) {
        return(max(v) + log(sum(exp(v - max(v)))))
    } else {
        return(-Inf)
    }
}

logsum1pexp <- function(w) {
    logic <- w > 0
    result <- vector("numeric", length(w))
    result[logic] <- w[logic] + log(1 + exp(-w[logic]))
    result[!logic] <- log(1 + exp(w[!logic]))
    return(result)
    # w may be a vector, but unlike logsumexp (where we get scalar
    # log(sum(exp(v)))), we'll get a vector log(exp(0)+exp(w))=log(1+exp(w)) as
    # a result return(sapply(w, function(x) logsumexp(c(0, x))))
}
