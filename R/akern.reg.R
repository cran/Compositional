akern.reg <- function(xnew, y, x, a = seq(0.1, 1, by = 0.1), h = seq(0.1, 1, length = 10) ) {

  est <- list()
  if ( min(y) == 0 )  a <- a[a > 0]
  la <- length(a)
  nh <- length(h)
  if ( !is.matrix(xnew) )  xnew <- as.matrix(xnew)
  nu <- dim(xnew)[1]
  D <- dim(y)[2]
  names <- paste("alpha", a)
  est <- sapply(names, function(x) NULL)
  di <-  -0.5 * Rfast::dista( xnew, x, square = TRUE)

  for ( i in 1:la ) {
    if ( abs( a[i] ) < 1e-9 ) {
      ua <- Compositional::alef(y, 0)$aff
	  for ( j in 1:nh ) {
        w <- exp( di / h[j] )  
        es <- exp( w %*% ua ) 
        est[[ i ]][[ j ]] <- es / Rfast::rowsums(es)       
      }
	} else {
      ua <- y^a[i]  
      ua <- ua / Rfast::rowsums(ua)
      for ( j in 1:nh ) {
        w <- exp( di / h[j] )  
        es <- ( w %*% ua )^( 1 / a[i] ) 
        est[[ i ]][[ j ]] <- es / Rfast::rowsums(es)       
      }
	}  
  }  ##  end  for ( i in 1:la ) {

  est
}
