################################################################################
#
# File name: fmt.R
#
# Author: Henrik Imberg
#
# Last edited: 2023-03-17
#
# Description: Simple text formatting function used for Table 2–4 (Section 6.2–6.4).
#
# INPUT: numeric vector x. 
#
# OUTPUT: character vector y. 
#
################################################################################

fmt <- function(x) { 
  if ( is.na(x) ) {
    y <- ""
  } else if( x > (1 - 1e-9) ) { 
    y <- "\\textbf{1.00}"
  } else if ( x > 0.995 ) {
    y <- "$>$0.99"
  } else {
    y <- sprintf("%.2f", x)
  }
  return(y)
}

fmt <- Vectorize(fmt)