calculate_sampling_scheme <- function(init, 
                                      grads, 
                                      hess,
                                      w = rep(1, length(init)), 
                                      design = c("MULT", "PO-WR", "PO-WOR"), 
                                      crit = c("A", "D", "E", "L"), 
                                      L = NULL, 
                                      maxiter = 1e2, 
                                      maxinner = 0,
                                      reltol = 1e-3, 
                                      abstol = 1e-3, 
                                      ftol = 1e-15,
                                      verbose = TRUE, 
                                      plot = TRUE) {
  
  design <- match.arg(design)
  crit <- match.arg(crit)
  force(w)
  
  # Some constants. 
  I <- solve(hess)
  N <- length(init)
  n <- sum(init)
  conv <- FALSE
  i <- 0
  jcount <- 0
  w_interp <- 0.01^(1 / maxinner) 
  
  # Init.
  mu <- init
  Gamma <- acov(mu, grads, hess, w, design)
  L <- chol_deriv(Gamma, crit, L)
  ci <- vapply(1:N, function(ix) sqrt(crossprod(L %*% I %*% grads[, ix])), numeric(1))
  val <- Phi(Gamma, crit, L)
  
  # Check convergence.
  if ( design == "PO-WOR" ) {
    ix <- which(mu >= 1) 
    jx <- which(mu < 1) 
    if ( length(ix) > 0 ) {
      cond1 <- abs(min(log(ci[ix]) - log(mu[ix])) - max(log(ci[jx]) - log(mu[jx])))
    } else {
      cond1 <- 1
    } 
    cond2 <- abs(max(log(ci[jx]) - log(mu[jx])) - min(log(ci[jx]) - log(mu[jx])))
    cond3 <- abs(max(log(ci[jx]) - log(mu[jx])) / min(log(ci[jx]) - log(mu[jx])))
    print(sprintf("cond1: %.6f, cond2: %.6f, cond3: %.6f", cond1, cond2, cond3))
    conv <- (cond1 < abstol) & ((cond2 < abstol) | (cond3 < reltol))
  } else {
    conv <- abs(max(log(ci) - log(mu)) - min(log(ci) - log(mu))) < abstol
    print(sprintf("Slack: %.2f", abs(max(log(ci) - log(mu)) - min(log(ci) - log(mu)))))
  }
  
  # Print.
  if ( verbose ) { print(sprintf("Iteration %d, objective function value: %.6f", 0, val)) }
  
  # Outer iteration.
  while ( !conv & i <= maxiter ) {
    i <- i + 1
    
    # Inner iteration.
    for ( j in 1:(maxinner + 1) ) {
      
      # Update sampling scheme.
      # Use grid search if no improvement was made in first inner iteration.
      if ( j == 1 ) {
        mu_tmp <- n * ci / sum(ci) 
      } else {
        mu_tmp <- (1 - w_interp) * mu + w_interp * mu_tmp 
        jcount <- jcount + 1
      }
      if ( design == "PO-WOR" ) {
        mu_next <- inclusionprobabilities(mu_tmp, n)
      } else {
        mu_next <- mu_tmp
      }
      
      # Update objective function etc.
      Gamma <- acov(mu_next, grads, hess, w, design)
      L <- chol_deriv(Gamma, crit, L)
      val_next <- Phi(Gamma, crit, L) 
      ci <- vapply(1:N, function(ix) sqrt(crossprod(L %*% I %*% grads[, ix])), numeric(1))
      
      # Break inner iteration if objective function value has improved.
      if ( val_next < val ) { break }
      
    }
    
    # If algorithm failed to improve objective function value.
    # if ( val_next > val ) {
    #   mu_next <- mu
    #   val_next <- val
    # }
    
    if ( verbose ) { print(sprintf("Iteration %d, objective function value: %.6f", i, val_next)) }
    
    # Check convergence.
    if ( design == "PO-WOR" ) {
      ix <- which(mu_next >= 1) 
      jx <- which(mu_next < 1) 
      if ( length(ix) > 0 ) {
        cond1 <- abs(min(log(ci[ix]) - log(mu_next[ix])) - max(log(ci[jx]) - log(mu_next[jx])))
      } else {
        cond1 <- 1
      } 
      cond2 <- abs(max(log(ci[jx]) - log(mu_next[jx])) - min(log(ci[jx]) - log(mu_next[jx])))
      cond3 <- abs(max(log(ci[jx]) - log(mu_next[jx])) / min(log(ci[jx]) - log(mu_next[jx])))
      print(sprintf("cond1: %.6f, cond2: %.6f, cond3: %.6f", cond1, cond2, cond3))
      conv <- (cond1 < abstol) & ((cond2 < abstol) | (cond3 < reltol))
    } else {
      conv <- abs(max(log(ci) - log(mu_next)) - min(log(ci) - log(mu_next))) < abstol
      print(sprintf("Slack: %.2f", abs(max(log(ci) - log(mu_next)) - min(log(ci) - log(mu_next)))))
    }
    
    # Stop for convergence.
    if ( conv | abs((val - val_next) / val) < ftol ) {
      mu <- mu_next
      val <- val_next
      break 
    }
    
    mu <- mu_next
    val <- val_next
    
  }
  
  if ( conv ) { flag <- 0 }
  else if ( val_next >  val | abs((val - val_next) / val) < ftol ) { flag <- 1 }
  else if ( i == maxiter ) { flag <- 2 }
  
  if ( flag == 0 ) {
    if ( i == 0 ) {
      message <- "Initial point is local minimum."
    } else {
      message <- "Algorithm converged."
    }
  }
  if ( flag == 1 ) {
    message <- "Warning: Failed to improve objective function value."
  }
  if ( flag == 2 ) {
    message <- "Did not converge. Maximal number of iterations reached."
  }
  if ( verbose ) { print(message) }
  
  if ( plot ) {
    if ( design == "PO-WOR" ) {
      jx <- which(mu < 1) 
    } else {
      jx <- 1:N
    }
    # slope <- max(mu[jx]) / max(ci[jx])
    plot(log(ci), log(mu), bty = "l")
    # abline(0, slope)
  }
  
  return(list(mu = mu, ci = ci, val = val, niter = i, ninner = jcount, status = flag, message = message))
  
}
