chol_deriv <- function(CovMat, 
                       crit = c("A", "D", "E", "L"), 
                       L = NULL) {
  
  crit <- match.arg(crit)
  
  if ( crit == "A" ) {
    L <- diag(1, nrow(CovMat))
  } else if ( crit == "D" ) {
    L <- sqrtm(solve(CovMat)) # Matrix square root.
  } else if ( crit == "E" ) {
    L <- eigen(CovMat)$vectors[, 1]
  } 
  
  return(L)
  
}