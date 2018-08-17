### test on convergence of algorithm

converged <- function(m_old, m_new, predicterror_old, predicterror_new, current.it, tol, iter.max) {
  
  maxchange <- max(0, abs(m_old-m_new), abs(predicterror_old - predicterror_new))
  
  if (!is.na(maxchange)) {
    if (current.it == iter.max) {
      # cat("reached max iterations (", current.it, ") and max change of parameters: ", round(maxchange, 3), "\n", sep="")
      # cat("iterations maxed \n")
    }
    if (maxchange < tol) {
      # cat("converged!\n")
    }
    return(maxchange < tol || current.it == iter.max)
  } else { # maxchange is NA
    cat("EP crashed!\n")
    return(TRUE)
  }
  
}