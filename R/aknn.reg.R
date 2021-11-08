## Response is compositional data
aknn.reg <- function(xnew, y, x, a = seq(0.1, 1, by = 0.1), k = 2:10, apostasi = "euclidean", rann = FALSE) {

  est <- list()
  if ( min(y) == 0 )  a <- a[a > 0]
  la <- length(a)
  nk <- length(k)
  if ( !is.matrix(xnew) )  xnew <- as.matrix(xnew)
  nu <- dim(xnew)[1]
  D <- dim(y)[2]
  names <- paste("alpha", a)
  est <- sapply(names, function(x) NULL)
  
  if ( rann ) {
    di <- RANN::nn2( data = x, query = xnew, k = max(k) )$nn.idx
  } else  di <- Rfast::dista( xnew, x, type = apostasi, k = max(k), index = TRUE, square = TRUE )
  
  for ( i in 1:la ) est[[ i ]] <- Compositional::frechet2( y, di, a[i], k )
    
  est
}
