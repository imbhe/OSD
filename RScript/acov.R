################################################################################
#
# File name: acov
#
# Author: Henrik Imberg
#
# Last edited: 2023-04-05
#
# Description: Evaluate approximate covariance matrix (Eq. (7)-(9), Section 2.2).
#
# INPUT: 
#   - mu:     sampling scheme
#   - grads:  gradients scheme
#   - hess:   Hessian
#   - design: sampling design
#
# OUTPUT: approximate covariance matrix.
#
################################################################################

acov <- function(mu, 
                 grads, 
                 hess, 
                 design = c("MULTI", "PO-WR", "PO-WOR")) {
  
  design <- match.arg(design)
  
  if ( any(mu < 0) ) { stop("mu < 0 not allowed.") }
  if ( design == "PO-WOR" & any(mu > 1) ) { stop("mu > 1 not allowed when design = PO-WOR.") }
  if ( design == "MULTI" & abs(round(sum(mu)) - sum(mu)) > 1e-8 ) { stop("Sample size must be integer when design = MULTI.") }
  
  N <- length(mu)
  n <- sum(mu)
  d <- nrow(grads)
  I <- solve(hess)
  cov <- matrix(0, d, d)
  
  if ( design %in% c("MULTI", "PO-WR") ) {
    for ( ix in 1:N ) {
      cov <- cov + tcrossprod(grads[, ix]) / mu[ix]
     }
  } else if ( design == "PO-WOR" ) {
    for ( ix in 1:N ) {
      cov <- cov +  tcrossprod(grads[, ix]) * (1 - mu[ix]) / mu[ix]
    }
  } 

  return(I %*% cov %*% I)
  
}