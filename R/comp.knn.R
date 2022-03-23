comp.knn <- function(xnew, x, ina, a = 1, k = 5, apostasi = "ESOV", mesos = TRUE) {
  n <- dim(x)[1]
  p <- dim(x)[2]
  ina <- as.numeric(ina)
  xnew <- as.matrix(xnew)
  xnew <- matrix(xnew, ncol = p ) ## makes sure xnew is a matrix
  nc <- max(ina)  ## The number of groups
  nu <- dim(xnew)[1]

  if ( !is.null(a) ) {
    if (apostasi == "CS"  &  a == 0 )  apostasi = "Ait"
  }

  if ( (apostasi == "taxicab" | apostasi == "Ait" | apostasi == "Hellinger") ) {

    if ( apostasi == "taxicab" ) {
	  if ( is.null(a) ) {
        zx <- x
        znew <- xnew
      } else {
        xa <- x^a
        zx <- xa / Rfast::rowsums( xa )  ## The power transformation is applied
        za <- xnew^a
        znew <- za / Rfast::rowsums( za )  ## The power transformation is applied
      }
	    g <- Rfast::knn(znew, ina, zx, k = k, dist.type = "manhattan", type = "C", freq.option = 1)
    } else if ( apostasi == "Ait" ) {
      xa <- Rfast::Log(x)
      zx <- xa - Rfast::rowmeans( xa )
      za <- Rfast::Log(xnew)
      znew <- za - Rfast::rowmeans( za )
	    g <- Rfast::knn(znew, ina, zx, k = k, dist.type = "euclidean", type = "C", freq.option = 1)

    } else if ( apostasi == "Hellinger" ) {
      g <- Rfast::knn(sqrt(xnew), ina, sqrt(x), k = k, dist.type = "euclidean", type = "C", freq.option = 1)
    }

  } else {
    ## all other methods
    klen <- length(k)
    g <- matrix(nrow = nu, ncol = klen)
    colnames(g) <- paste("k=", k, sep = "")

    if ( apostasi == "angular" ) {
      zx <- sqrt(x)
      znew <- sqrt(xnew)
      disa <- tcrossprod(zx, znew )
      disa[disa >= 1] <- 1
      disa[ disa <=  -1 ] <-  -1
      disa <- acos(disa)

    } else if ( apostasi == "ESOV" ) {
	    if ( is.null(a) ) {
        zx <- x
        znew <- xnew
	    } else {
        xa <- x^a
        zx <- xa / Rfast::rowsums( xa )  ## The power transformation is applied
        za <- xnew^a
        znew <- za / Rfast::rowsums( za )  ## The power transformation is applied
	    }
      tzx <- t(zx)
      disa <- matrix( nrow = nu, ncol = dim(tzx)[2] )
      for (i in 1:nu) {
        zan <- znew[i, ]
        ma <- 0.5 * ( tzx + zan )
        disa[i, ] <- colSums( zan * log( zan / ma ) + tzx * log( tzx/ma ), na.rm = TRUE )
      }

    } else if ( apostasi == "CS" ) {
      if ( is.null(a) ) {
        zx <- x
        znew <- xnew
      } else {
	    xa <- x^a
        zx <- xa / Rfast::rowsums( xa )  ## The power transformation is applied
        za <- xnew^a
        znew <- za / Rfast::rowsums( za )  ## The power transformation is applied
      }
     	tzx <- t(zx)
     	disa <- matrix( nrow = nu, ncol = dim(tzx)[2] )
     	for (i in 1:nu) {
        znewi <- znew[i, ]
        sa <- ( tzx - znewi )^2 / ( tzx + znewi )
        sa[is.infinite(sa)] <- 0
        disa[i, ] <- Rfast::colsums( sa )
      }  ## end for (i in 1:nu)
        ## disa <- sqrt(disa) / abs(a) * sqrt(2 * p) not necessary to take the sqrt and then divide and multiply with constants everywhere
    }  ## end if (apostasi == "CS")

    disa <- Rfast::rowOrder(disa)
    for (j in 1:klen) {
       for (i in 1:nu) {
         mod <- table( ina[ which(disa[i, ] <= k[j]) ] )
         g[i, j] <- as.numeric( names(mod)[ which.max(mod) ] )
       }
    }

  }  ## end of other methods

  g
}
