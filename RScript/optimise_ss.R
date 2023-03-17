################################################################################
#
# File name: optimise_ss.R
#
# Author: Henrik Imberg
#
# Last edited: 2023-03-17
#
# Description: Find optimal sampling scheme using the iterative L-approximation 
#              algorithm (Algorithm 2, Section 3.3).
#
# INPUT: 
#   - n           subsample size.
#   - grads       gradients scheme.
#   - hess        Hessian.
#   - design      sampling design.
#   - crit        optimality criterion.
#   - L           L-matrix for linear optimality criteria, including c-optimality. 
#   - r           r parameter for the Phi_r-optimality criterion.
#   - init        initial sampling scheme (defaults to uniform).
#   - maxiter     maximal number of iterations. 
#   - Delta       tolerance paramter for stationarity conditions.
#   - eps         tolerance paramter for relative improvement of objective function. 
#   - print_level 0 = no printed output, 1 = print after final iteration, 
#                 2 = print at each iteration, 3 = detailed output at each iteration.
#   - plot        Should plots be produced? 
#   - ggplot      Should ggplot be produced?
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

optimise_ss <- function(n,
                 grads, 
                 hess,
                 design = c("MULTI", "PO-WR", "PO-WOR"), 
                 crit = c("A", "c", "D", "E", "L", "Phi_r"), 
                 L = NULL, 
                 r = NULL,
                 init = rep(1, ncol(grads)), 
                 maxiter = 1e2, 
                 Delta = 1e-1, 
                 eps = 1e-3, 
                 print_level = 1, 
                 plot = FALSE, 
                 ggplot = FALSE) {
  
  library("expm")
  library("sampling")
  source("acov.R")
  source("Phi.R")
  
  if ( any(init < 0) ) { stop("Initial value < 0 not allowed.") }
  if ( length(init) != ncol(grads) ) { stop("Dimensions do not agree. length(init) != ncol(Y).") }
  if ( !is.na(r) && !is.null(r) && r <= 0 ) { stop("r must be > 0.") }
  if ( ggplot == TRUE ) { plot = TRUE }
  
  if(print_level > 0) {
    print_crit <- ifelse(crit == "Phi_r", paste0("Phi_", r), crit)
    cat(sprintf("--- %s-optimality criterion ---\n", print_crit))
  }
  
  t0 <- Sys.time()
  design <- match.arg(design)
  crit <- match.arg(crit)

  # Initialisation.
  ci <- init^2
  I <- solve(hess) # Inverse Hessian.
  N <- length(init)
  vals <- rep(NA, maxiter + 1)
  delta <- rep(NA, maxiter + 1)
  frel <- NULL

  # Iterate.
  for ( i in 1:(maxiter + 1) ) {
    
    # Update sampling scheme.
    mu <- n * sqrt(ci) / sum(sqrt(ci)) 
    if ( design == "PO-WOR" ) {
      mu <- inclusionprobabilities(mu, n) # Inclusion probabilities in (0, 1]
    }

    # Evaluate objective function etc. 
    Gamma <- acov(mu, grads, hess, w, design)
    val <- Phi(Gamma, crit, L, r)
    vals[i] <- val

    # Evaluate L-matrix.
    if ( crit == "A" ) {
      L <- diag(1, nrow(Gamma))
    } else if ( crit == "D" ) {
      L <- sqrtm(solve(Gamma)) # Matrix square root.
    } else if ( crit == "E" ) {
      L <- eigen(Gamma)$vectors[, 1]
    } else if ( crit == "Phi_r" ) {
      eig <- eigen(Gamma)
      P <- eig$vectors
      D <- diag(eig$values^(0.5 * (r - 1)))
      L <- P %*% D %*% t(P)
    }
    
    # Calculate ci's.
    LtI <- t(L) %*% I
    ci <- vapply(1:N, function(ix) crossprod(LtI %*% grads[, ix]), numeric(1))
  
    # Check stationarity conditions.
    if ( design == "PO-WOR" ) {
      ix <- which(mu >= 1) # E = {i in D: mui >= 1}.
      jx <- which(mu < 1) # D \ E.
      delta1 <- abs(max(0.5 * log(ci[jx]) - log(mu[jx])) - min(0.5 * log(ci[jx]) - log(mu[jx])))
      if ( length(ix) > 0 ) {
        delta2 <- max(0.5 * log(ci[jx]) - log(mu[jx])) - min(0.5 * log(ci[ix]) - log(mu[ix]))
      } else {
        delta2 <- 0
      } 
      delta[i] <- max(c(delta1, delta2))
    } else {
      delta[i] <- abs(max(0.5 * log(ci) - log(mu)) - min(0.5 * log(ci) - log(mu)))
    }
    
    # Improvement in objective function.
    if ( i > 0 ) {
      frel <- abs(vals[i] - vals[i - 1]) / vals[i - 1]
    }

    # Print.
    if ( print_level > 1 ) { 
      cat(sprintf("Iteration %d, objective = %g.\n", i - 1, val))
    } 
    if ( print_level > 2 ) {
      cat(sprintf("Stationarity conditions: %g < %g is %s.\nRelative improvement of objective function: %g < %g is %s.\n", 
                  delta[i], Delta, ifelse(delta[i] < Delta, "TRUE", "FALSE"),
                  frel, eps, ifelse(frel < eps, "TRUE", "FALSE")))
    }

    # Stop if algorithm has converged/diverged.
    iters <- ifelse(i == 2, "iteration", "iterations")
    if ( delta[i] < Delta ) {
      flag <- 0
      if ( i == 1) {
        message <- sprintf("Algorithm converged in %d %s. Intial sampling scheme is a stationary point.", i - 1, iters)
      } else {
        message <- sprintf("Algorithm converged in %d %s. Stationarity conditions fulfilled.", i - 1, iters)
      }
      break
    } else if ( i > 1 && frel < eps ) {
      flag <- 0
      message <- sprintf("Algorithm converged in %d %s. Relative improvement in objective function less than %g.", i - 1, iters, eps)
      break
    } else if ( delta[i] != min(delta, na.rm = TRUE) ) {
      flag <- 1
      message <- sprintf("Algorithm diverged after %d %s. Stationarity conditions worsened during iteration.", i - 1, iters)
      break
    } else if ( val != min(vals, na.rm = TRUE) ) {
      flag <- 1
      message <- sprintf("Algorithm diverged after %d %s. Objective function increased during iteration.", i - 1, iters)
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
    cat(sprintf(" %.2f sec elapsed.\n\n", as.numeric(t1 - t0))) 
  }
  
  # Plot.
  if ( plot & flag != 1 ) {
    if ( design == "PO-WOR" ) {
      jx <- which(mu < 1) 
    } else {
      jx <- 1:N
    }
    
    if ( ggplot ) {
      
      library("cowplot")
      library("ggplot2")
      
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
              niter = i - 1, 
              status = flag, 
              message = message,
              Gamma = Gamma,
              elapsed = as.numeric(t1 - t0)))
  
}