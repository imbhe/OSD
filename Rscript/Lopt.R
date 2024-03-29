################################################################################
#
# File name: lopt.R
#
# Author: Henrik Imberg
#
# Last edited: 2023-04-05
#
# Description:  Find L-optimal sampling scheme using Algorithm 1, Section 3.3.
#
# INPUT: 
#   - n           subsample size.
#   - grads       gradients scheme.
#   - hess        Hessian.
#   - design      sampling design.
#   - L           L-matrix for linear optimality criteria, 
#                 including A- and c-optimality. 
#
# OUTPUT: 
#   - mu          Sampling scheme at final iteration.
#   - elapsed     Time elapsed. 
#
################################################################################

lopt <- function(n,
                 grads, 
                 hess,
                 L = NULL,
                 design = c("MULTI", "PO-WR", "PO-WOR"), 
                 print_level = 1) {
  
  t0 <- Sys.time()
  
  library("sampling")
  
  if ( is.null(L) ) { stop('Real matrix L must be specified.') }
  
  if(print_level > 0) {
    cat(sprintf("--- L-optimality criterion ---\n"))
  }
  
  design <- match.arg(design)
  
  # Calculate ci's.
  I <- solve(hess) # Inverse Hessian.
  LtI <- t(L) %*% I
  ci <- vapply(1:N, function(ix) crossprod(LtI %*% grads[, ix]), numeric(1))
  
  # Calculate sampling scheme.
  mu <- n * sqrt(ci) / sum(sqrt(ci)) 
  if ( design == "PO-WOR" ) {
    mu <- inclusionprobabilities(mu, n) # Inclusion probabilities in (0, 1]
  }
  
  t1 <- Sys.time()
  
  # Print.
  if ( print_level > 0 ) { 
    cat(sprintf(" %.2f sec elapsed.\n\n", as.numeric(t1 - t0))) 
  }
  
  return(list(mu = mu, 
              niter = NA, 
              elapsed = as.numeric(t1 - t0)))
  
}