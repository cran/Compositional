## computes the ESOV distance between two vectors, x1 and x2
es <- function(x1, x2)  sum( x1 * log( 2 * x1 / (x1 + x2) ) + x2 * log( 2 * x2 / (x1 + x2) ), na.rm = TRUE )
