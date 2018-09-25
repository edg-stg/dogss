my_cvSGL <- function(data, index, nfolds=10, standardize=FALSE, thresh=0.001, lambdas=NULL, maxit=1000) { # 
  Y <- data$y
  X <- data$x
  nobs <- length(Y)
  foldids <- sample(rep(seq(nfolds), length = nobs))
  outlist <- as.list(seq(nfolds))
  
  main_run <- SGL(data=data, index=index, standardize=standardize, thresh = thresh, lambdas=lambdas, maxit=maxit) # 
  
  lambdas <- main_run$lambdas
  
  predmat <- matrix(NA, length(Y), length(lambdas))
  betas <- vector("list", nfolds) # matrix(NA, dim(X)[2], nfolds)

  for (i in seq(nfolds)) {
    which <- foldids == i
    Y_sub <- Y[!which]
    X_sub <- X[!which, , drop=FALSE]
    outlist[[i]] = SGL(data=list(x=X_sub, y=Y_sub), index=index, standardize=standardize, lambdas=lambdas, thresh = thresh, maxit=maxit) # 
    betas[[i]] <- outlist[[i]]$beta
    preds <- matrix(NA, nrow=sum(which), ncol=length(lambdas))
    for (l in seq(lambdas)) {
      preds[, l] <- X[which, , drop = FALSE] %*% betas[[i]][, l]
    }
    predmat[which, ] <- preds
  }

  cvraw <- (Y - predmat)^2
  cvm <- apply(cvraw, 2, mean, na.rm = TRUE)
  cvsd <- sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, mean, na.rm = TRUE)/(nobs - 1))

  which_cvm_min <- which.min(cvm)
  which_1se <- min(which(cvm<=min(cvm)+cvsd[which_cvm_min]))

  beta_final <- main_run$beta[, which_1se]

  out <- list(beta=beta_final, lambda_1se=lambdas[which_1se], fit=main_run) # cvm = cvm, cvsd = cvsd,

  return(out)
}
