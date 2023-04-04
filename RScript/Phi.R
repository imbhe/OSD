################################################################################
#
# File name: Phi.R
#
# Author: Henrik Imberg
#
# Last edited: 2023-03-17
#
# Description: Evaluate Phi-optimality objective function (Table 1, Section 2.3). 
#
# INPUT: 
#   - CovMat      (Approximate) covariance matrix.
#   - crit        optimality criterion.
#   - L           L-matrix for linear optimality criteria, including c-optimality. 
#   - q           q parameter for the Phi_q-optimality criterion.
#
# OUTPUT: value of objective function Phi.
#
################################################################################

Phi <- function(CovMat, 
                crit = c("A", "c", "D", "E", "L", "Phi_q"), 
                L = NULL, 
                q = NULL) {
  
  crit <- match.arg(crit)
  if ( !is.na(q) && !is.null(q) && q <= 0 ) { stop("q must be > 0.") }
  if ( crit == "Phi_q" && is.null(q) ) { stop('Real number q > 0 must be specified when crit = "Phi_q".') }
  if ( crit == "L" && is.null(L) ) { stop('Real matrix L must be specified when crit = "L".') }

  d <- ncol(CovMat)
  
  if ( crit == "A" ) {
    val <- sum(diag(CovMat))
  } else if ( crit == "D" ) {
    val <- as.numeric(exp(determinant(CovMat, log = TRUE)$modulus / ncol(CovMat)))
  } else if ( crit == "E" ) {
    val <- max(eigen(CovMat)$values)
  } else if ( crit %in% c("c", "L") ) {
    val <- sum(diag(CovMat %*% tcrossprod(L)))
  } else if ( crit == "Phi_q" ) {
    eig <- eigen(CovMat)
    P <- eig$vectors
    D <- diag(eig$values^q)
    val <- sum(diag(P %*% D %*% t(P)))^{1/q}
  }
  
  return(val)
  
}