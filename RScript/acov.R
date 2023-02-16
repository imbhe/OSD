acov <- function(mu, 
                 grads, 
                 hess, 
                 w = rep(1, length(mu)), 
                 design = c("MULT", "PO-WR", "PO-WOR")) {
  
  design <- match.arg(design)
  
  if ( any(mu < 0) ) { stop("mu < 0 not allowed.") }
  if ( design == "PO-WOR" & any(mu > 1) ) { stop("mu > 1 not allowed when design = PO-WOR.") }
  if ( design == "MULT" & abs(round(sum(mu)) - sum(mu)) > 1e-8 ) { stop("Sample size must be integer when design = MULT.") }
  
  N <- length(mu)
  n <- sum(mu)
  d <- nrow(grads)
  I <- solve(hess)
  cov <- matrix(0, d, d)
  
  if ( design == "MULT" ) {
    for ( ix in 1:N ) {
      cov <- cov + tcrossprod(grads[, ix]) / mu[ix]
      # Note: Formula simplifies since increment from loop below cancels out.
      # for ( jx in 1:N ) {
      #   cov <- cov - tcrossprod(grads[, ix], grads[, jx]) / n
      # }
    }
  } else if ( design == "PO-WOR" ) {
    for ( ix in 1:N ) {
      cov <- cov +  tcrossprod(grads[, ix]) * (1 - mu[ix]) / mu[ix]
    }
  } else if ( design == "PO-WR") {
    for ( ix in 1:N ) {
      cov <- cov + tcrossprod(grads[, ix]) / mu[ix]
    }
  } 

  return(I %*% cov %*% I)
  
}