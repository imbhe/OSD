Phi <- function(CovMat, 
                crit = c("A", "D", "E", "L", "Phi_p"), 
                L = NULL, 
                p = NULL) {
  
  crit <- match.arg(crit)
  if ( !is.null(p) && p <= 0 ) { stop("p must be > 0.") }
  
  if ( crit == "A" ) {
    val <- sum(diag(CovMat))
  } else if ( crit == "D" ) {
    val <- log(det(CovMat))
  } else if ( crit == "E" ) {
    val <- max(eigen(CovMat)$values)
  } else if ( crit == "L" ) {
    val <-sum(diag(CovMat %*% tcrossprod(L)))
  } else if ( crit == "Phi_p" ) {
    eig <- eigen(CovMat)
    P <- eig$vectors
    D <- diag(eig$values^p)
    val <- sum(diag(P %*% D %*% t(P)))^{1/p}
  }
  
  return(val)
  
}