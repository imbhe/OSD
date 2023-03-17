opt_sampling_scheme <- function(n,
                                grads, 
                                hess,
                                w = rep(1, length(init)), 
                                design = c("MULT", "PO-WR", "PO-WOR"), 
                                crit = c("A", "c", "D", "E", "L", "Phi_p"), 
                                L = NULL, 
                                p = NULL,
                                init = rep(1, ncol(grads)), 
                                maxiter = 50, 
                                tol = 1e-1, 
                                freltol = 1e-3, 
                                kappa_tol = 1e9,
                                print_level = 0, 
                                plot = FALSE, 
                                ggplot = FALSE) {
  
  require("expm")
  
  if ( any(init < 0) ) { stop("Initial value < 0 not allowed.") }
  if ( length(init) != ncol(grads) ) { stop("Dimensions do not agree. length(init) != ncol(Y).") }
  if ( !is.na(p) && !is.null(p) && p <= 0 ) { stop("p must be > 0.") }
  
  tic()
  design <- match.arg(design)
  crit <- match.arg(crit)
  
  # Initialisation.
  I <- solve(hess)
  N <- length(init)
  conv <- FALSE
  i <- 0
  vals <- rep(NA, maxiter + 1)
  
  # Initial sampling scheme of (expected) size n.
  mu <- n * init / sum(init) 
  if ( design == "PO-WOR" & any(mu > 1) ) { # Make sure initial sampling scheme is feasible for PO-WOR.
    mu <- inclusionprobabilities(mu, n) 
  }
  
  # Iterate.
  while ( 1 ) {
    
    # Evaluate objective function etc.
    Gamma <- acov(mu, grads, hess, w, design)
    val <- Phi(Gamma, crit, L, p)
    vals[i + 1] <- val
    if ( kappa(Gamma) > kappa_tol ) { # Stop if problem too ill-conditioned.
      break
    }
    
    # Evaluate L-matrix.
    if ( crit == "A" ) {
      L <- diag(1, nrow(Gamma))
    } else if ( crit == "D" ) {
      L <- sqrtm(solve(Gamma)) # Matrix square root.
    } else if ( crit == "E" ) {
      L <- eigen(Gamma)$vectors[, 1]
    } else if ( crit == "Phi_p" ) {
      eig <- eigen(Gamma)
      P <- eig$vectors
      D <- diag(eig$values^(0.5 * (p - 1)))
      L <- P %*% D %*% t(P)
    }
    LtI <- t(L) %*% I
    ci <- vapply(1:N, function(ix) sqrt(crossprod(LtI %*% grads[, ix])), numeric(1))
    
    # Check convergence.
    if ( design == "PO-WOR" ) {
      ix <- which(mu >= 1) 
      jx <- which(mu < 1) 
      cond1 <- abs(max(log(ci[jx]) - log(mu[jx])) - min(log(ci[jx]) - log(mu[jx])))
      if ( length(ix) > 0 ) {
        cond2 <- max(log(ci[jx]) - log(mu[jx])) - min(log(ci[ix]) - log(mu[ix]))
      } else {
        cond2 <- 0
      } 
      conv <- (cond1 < tol) & (cond2 < tol)
      if ( print_level > 1 ) {
        note <- sprintf("cond1 %.6f < %.3f is %s. cond2 %.6f < %.3f is %s. kappa = %.f.", 
                        cond1, tol, ifelse(cond1 < tol, "TRUE", "FALSE"),
                        cond2, tol, ifelse(cond2 < tol, "TRUE", "FALSE"), 
                        kappa(Gamma))
      }
    } else {
      cond1 <- abs(max(log(ci) - log(mu)) - min(log(ci) - log(mu)))
      conv <- cond1 < tol  
      if ( print_level > 1 ) {
        note <- sprintf("cond %.6f < %.3f is %s. kappa = %.f.", 
                        cond1, tol, ifelse(cond1 < tol, "TRUE", "FALSE"),
                        kappa(Gamma))
      }
    }
    
    # Print.
    if ( print_level == 1 ) { 
      cat(sprintf("Iteration %d. objective = %.6f.\n", i, val))
    } else if ( print_level == 2 ) { 
      cat(sprintf("Iteration %d. objective = %.6f. %s\n", i, val, note))
    }
    
    # Stop if algorithm has converged.
    if ( conv | i == maxiter | maxiter == 0 ) { break }
    
    # Else: update sampling scheme.
    i <- i + 1    
    mu <- n * ci / sum(ci) 
    if ( design == "PO-WOR" ) {
      mu <- inclusionprobabilities(mu, n)
    }
    
  }
  
  t <- toc(quiet = TRUE)
  
  # Flags.
  if ( conv ) { 
    flag <- 0 
  }  else if ( kappa(Gamma) > kappa_tol ) {
    flag <- 1
  } else { 
    flag <- 2 
  }
  
  # Messages.
  if ( flag == 0 ) {
    minval <- min(vals[1:i])
    if ( (val - minval) / abs(minval) > freltol ) {
      message <- "Algorithm converged but value of objective function has increased. Problem may be ill-conditioned."
    } else if ( i == 0 ) {
      message <- "Initial point is local minimum."
    } else {
      message <- "Algorithm converged."
    }
  } else if ( flag == 1 ) {
    message <- sprintf("Algorithm terminated. Problem too ill-conditioned (kappa = %.0f).", kappa(Gamma))
  } else if ( flag == 2 ) {
    message <- "Did not converge. Maximal number of iterations reached."
  } 
  
  # Print.
  if ( print_level > 0 ) { 
    cat(message) 
    cat(sprintf(" %.2f sec elapsed.\n", as.numeric(t$toc - t$tic))) 
  }
  
  # Plot.
  if ( plot & flag != 1 ) {
    if ( design == "PO-WOR" ) {
      jx <- which(mu < 1) 
    } else {
      jx <- 1:N
    }
    
    if ( ggplot ) {
      fig1 <- ggplot() + 
        geom_line(aes(x = 0:i, y = vals[1:(i + 1)])) + 
        geom_point(aes(x = 0:i, y = vals[1:(i + 1)])) + 
        labs(x = "Iteration",
             y = "Objective") + 
        theme_classic()
      
      fig2 <- ggplot() + 
        geom_abline(aes(intercept = log(max(mu[jx])) - log(max(ci[jx])), slope = 1)) + 
        geom_point(aes(x = log(ci), y = log(mu))) + 
        labs(x = expression(log(c[i])),
             y = expression(log(mu[i]))) + 
        theme_classic()
      
      print(plot_grid(fig1, fig2))
    } else {
      par(mfrow = c(1, 2))
      plot(0:i, vals[1:(i + 1)], type = "b", bty = "l", xlab = "Iteration", ylab = "Objective")
      plot(log(ci), log(mu), bty = "l", xlab = expression(log(c[i])), ylab = expression(log(mu[i])))
      abline(log(max(mu[jx])) - log(max(ci[jx])), 1)
      par(mfrow = c(1, 1))
    }
    
  }
  
  return(list(mu = mu, 
              val = val, 
              niter = i, 
              status = flag, 
              message = message,
              crit = crit, 
              design = design,
              L = L, 
              p = p,
              elapsed = as.numeric(t$elapsed)))
  
}
