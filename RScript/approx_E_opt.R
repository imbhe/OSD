approx_E_opt <- function(n,
                         grads, 
                         hess,
                         w = rep(1, length(init)), 
                         design = c("MULT", "PO-WR", "PO-WOR"), 
                         eff = 0.95,
                         init = rep(1, ncol(grads)), 
                         maxiter = 20,
                         tol = 1e-1, 
                         freltol = 1e-3, 
                         kappa_tol = 1e9,
                         print_level = 2, 
                         plot = TRUE, 
                         ggplot = FALSE) {
 
  d <- nrow(hess)
  p_target <- ceiling(-log(d) / log(eff))
  niter <- 0
  
  for ( p in 1:p_target ) {
    
    if ( p == 1) {
      new <- opt_sampling_scheme(n = n,
                                 grads = grads, 
                                 hess = hess,
                                 w = w,
                                 design = design,
                                 crit = "Phi_p",
                                 p = p,
                                 init = rep(n / N, N),
                                 maxiter = maxiter,
                                 tol = tol,
                                 freltol = freltol, 
                                 kappa_tol = kappa_tol,
                                 print_level = print_level, 
                                 plot = plot, 
                                 ggplot = ggplot)
    } else {
      new <- opt_sampling_scheme(n = n,
                                 grads = grads, 
                                 hess = hess,
                                 w = w,
                                 design = design,
                                 crit = "Phi_p",
                                 p = p,
                                 init = old$mu,
                                 maxiter = maxiter,
                                 tol = tol,
                                 freltol = freltol, 
                                 kappa_tol = kappa_tol,
                                 print_level = print_level, 
                                 plot = plot, 
                                 ggplot = ggplot)
    }
    
    
    niter <- niter + new$niter
    
    if ( new$status != 0 ) {
      break
    }
    
    old <- new
    
  }
  
  old$niter <- niter

  return(old)
  
}