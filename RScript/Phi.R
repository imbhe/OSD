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
#   - s           s parameter for the Phi_s-optimality criterion.
#
# OUTPUT: value of objective function Phi.
#
################################################################################

Phi <- function(CovMat, 
                crit = c("A", "c", "D", "E", "L", "Phi_s"), 
                L = NULL, 
                s = NULL) {
  
  crit <- match.arg(crit)
  if ( !is.na(s) && !is.null(s) && s <= 0 ) { stop("s must be > 0.") }
  if ( crit == "Phi_s" && is.null(s) ) { stop('Real number s > 0 must be specified when crit = "Phi_s".') }
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
  } else if ( crit == "Phi_s" ) {
    eig <- eigen(CovMat)
    P <- eig$vectors
    D <- diag(eig$values^s)
    val <- sum(diag(P %*% D %*% t(P)))^{1/s}
  }
  
  return(val)
  
}