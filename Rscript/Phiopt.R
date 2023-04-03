################################################################################
#
# File name: phiopt.R
#
# Author: Henrik Imberg
#
# Last edited: 2023-03-20
#
# Description:  Find Phi-optimal sampling scheme using the fixed point method (Algorithm 2, Section 3.3).
#
# INPUT: 
#   - n           subsample size.
#   - grads       gradients scheme.
#   - hess        Hessian.
#   - design      sampling design.
#   - crit        optimality criterion.
#   - L           L-matrix for linear optimality criteria, including c-optimality. 
#   - s           s parameter for the Phi_s-optimality criterion.
#   - init        initial sampling scheme (defaults to uniform).
#   - maxiter     maximal number of iterations. 
#   - eps         tolerance parameter for relative improvement of objective function. 
#   - print_level 0 = no printed output, 
#                 1 = print after final iteration, 
#                 2 = print objective function value at each iteration.
#
# OUTPUT: 
#   - mu          Sampling scheme at final iteration.
#   - val         Objective function value.
#   - niter       Number of iterations.
#   - status      0 = converged, 1 = diverged, 2 = maximal number of iteration reached.
#   - message     Status message.
#   - Gamma       Approximate covariance matrix at final iteration. 
#   - elapsed     Time elapsed. 
#
################################################################################

phiopt <- function(n,
                   grads, 
                   hess,
                   crit = c("A", "c", "D", "E", "L", "Phi_s"), 
                   L = NULL, 
                   s = NULL,
                   design = c("MULTI", "PO-WR", "PO-WOR"), 
                   init = rep(1, ncol(grads)), 
                   maxiter = 1e2, 
                   eps = 1e-3, 
                   print_level = 1) {
  
  t0 <- Sys.time()
  
  library("expm")
  library("sampling")
  source("Rscript/acov.R")
  source("Rscript/Phi.R")
  
  if ( any(init < 0) ) { stop("Initial value < 0 not allowed.") }
  if ( length(init) != ncol(grads) ) { stop("Dimensions do not agree. length(init) != ncol(Y).") }
  if ( !is.na(s) && !is.null(s) && s <= 0 ) { stop("s must be > 0.") }
  if ( crit == "Phi_s" && is.null(s) ) { stop('Real number s > 0 must be specified when crit = "Phi_s".') }
  if ( crit == "L" && is.null(L) ) { stop('Real matrix L must be specified when crit = "L".') }
  
  if(print_level > 0) {
    print_crit <- ifelse(crit == "Phi_s", paste0("Phi_", s), crit)
    cat(sprintf("--- %s-optimality criterion ---\n", print_crit))
  }
  
  design <- match.arg(design)
  crit <- match.arg(crit)
  
  # Initialisation.
  ci <- init^2
  I <- solve(hess) # Inverse Hessian.
  N <- length(init)
  vals <- rep(NA, maxiter + 1)
  frel <- NULL
  
  # Iterate.
  for ( i in 1:(maxiter + 1) ) {
    
    # Update sampling scheme.
    mu <- n * sqrt(ci) / sum(sqrt(ci)) 
    if ( design == "PO-WOR" ) {
      mu <- inclusionprobabilities(mu, n) # Inclusion probabilities in (0, 1]
    }
    
    # Evaluate objective function etc. 
    Gamma <- acov(mu, grads, hess, design)
    val <- Phi(Gamma, crit, L, s)
    vals[i] <- val
    
    # Evaluate L-matrix.
    if ( crit == "A" ) {
      L <- diag(1, nrow(Gamma))
    } else if ( crit == "D" ) {
      L <- sqrtm(solve(Gamma)) # Matrix square root.
    } else if ( crit == "E" ) {
      L <- eigen(Gamma)$vectors[, 1]
    } else if ( crit == "Phi_s" ) {
      eig <- eigen(Gamma)
      P <- eig$vectors
      D <- diag(eig$values^(0.5 * (s - 1)))
      L <- P %*% D %*% t(P)
    }
    
    # Calculate ci's.
    LtI <- t(L) %*% I
    ci <- vapply(1:N, function(ix) crossprod(LtI %*% grads[, ix]), numeric(1))
    
      # Improvement in objective function.
    if ( i > 0 ) {
      frel <- abs(vals[i] - vals[i - 1]) / vals[i - 1]
    }
    
    # Print.
    if ( print_level > 1 ) { 
      cat(sprintf("Iteration %d, objective = %g.\n", i - 1, val))
    } 

    # Stop if algorithm has converged/diverged.
    iters <- ifelse(i == 2, "iteration", "iterations")
    if ( i > 1 && frel < eps ) {
      flag <- 0
      message <- sprintf("Algorithm converged in %d %s", i - 1, iters)
      break
    } else if ( val != min(vals, na.rm = TRUE) ) {
      flag <- 1
      message <- sprintf("Algorithm diverged after %d %s. Value of objective function increased during iteration.", i - 1, iters)
      break
    } else if ( i == maxiter ) {
      flag <- 2
      message <- "Did not converge. Maximal number of iterations reached."
    }
    
  }
  
  t1 <- Sys.time()
  
  # Print.
  if ( print_level > 0 ) { 
    cat(message) 
    cat(sprintf(" %.2f sec elapsed.\n\n", as.numeric(difftime(t1, t0, units='secs'))))
  }

  return(list(mu = mu, 
              val = val, 
              niter = i - 1, 
              status = flag, 
              message = message,
              Gamma = Gamma,
              elapsed = as.numeric(difftime(t1, t0, units='secs'))))
  
}