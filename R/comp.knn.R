comp.knn <- function(xnew, x, ina, a = 1, k = 5, type = "S", apostasi = "ESOV", mesos = TRUE) {
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

  if ( (apostasi == "taxicab" | apostasi == "Ait" | apostasi == "Hellinger")  &  type == "S" ) {

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
	    g <- Rfast::knn(znew, ina, zx, k = k, dist.type = "mahattan", type = "C", freq.option = 1)
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
	    disa <- Rfast::dista(znew, zx, "manhattan", trans = FALSE)

    } else if ( apostasi == "Ait" ) {
      xa <- Rfast::Log(x)
      zx <- xa - Rfast::rowmeans( xa )
      za <- Rfast::Log(xnew)
      znew <- za - Rfast::rowmeans( za )
  	  disa <- Rfast::dista(znew, zx, trans = FALSE)

    } else if ( apostasi == "Hellinger" ) {
      disa <- Rfast::dista(sqrt(xnew), sqrt(x), "euclidean", trans = FALSE ,square = TRUE)

    } else if ( apostasi == "angular" ) {
      zx <- sqrt(x)
      znew <- sqrt(xnew)
      disa <- tcrossprod(zx, znew )
      disa[disa >= 1] <- 1
      disa[ disa <=  -1 ] <-  -1
      disa <- acos(disa)

    } else if ( apostasi == "ESOV" ) {
      g <- matrix(0, nu, klen)
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
      if (type == "NS") {
        disa <- matrix(0, n, nu)
        for (i in 1:nu) {
          zan <- znew[i, ]
          ma <- 0.5 * ( tzx + zan )
          disa[, i] <- colSums( zan * log( zan / ma ) + tzx * log( tzx / ma ), na.rm = TRUE )
        }
        ta <- matrix(nrow = nu, ncol = nc)
	  if (mesos) {
          for (j in 1:klen) {
            for (m in 1:nc) {
              apo <- disa[ina == m, ]
              apo <- Rfast::colSort(apo)
              ta[, m] <- Rfast::colmeans( apo[1:k[j], , drop = FALSE] )
            }
            g[, j] <- Rfast::rowMins(ta)
          }
	  } else {
	    for (j in 1:klen) {
            for (m in 1:nc) {
              apo <- disa[ina == m, ]
              apo <- Rfast::colSort(apo)
              ta[, m] <- Rfast::colhameans( apo[1:k[j], , drop = FALSE] )
            }
            g[, j] <- Rfast::rowMins(ta)
          }
	  } ## end if (mesos)
      } else {
        for (i in 1:nu) {
          zan <- znew[i, ]
          ma <- 0.5 * ( tzx + zan )
          di <- colSums( zan * log( zan / ma ) + tzx * log( tzx/ma ), na.rm = TRUE )
          di <- Rfast::Order(di)
          di <- di[ di <= max(k) ]
          for (j in 1:klen) {
            ind <- di[ 1:k[j] ]
            a <- Rfast::Table( ina[ind] )
            b <- as.numeric( names(a) )
            g[i, j] <- b[which.max(a)]
          }  ## end for (j in 1:klen)
        }  ## end for (i in 1:nu)
      }
    } else if ( apostasi == "CS" ) {
      g <- matrix(0, nu, klen)
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
      if (type == "NS") {
        disa <- matrix(0, n, nu)
        for (i in 1:nu) {
          znewi <- znew[i, ]
          sa <- ( tzx - znewi )^2 / ( zx + znewi )
          sa[is.infinite(sa)] <- 0
          disa[, i] <- Rfast::colsums( sa )
        }
        ta <- matrix(nrow = nu, ncol = nc)
	      if (mesos) {
          for (j in 1:klen) {
            for (m in 1:nc) {
              apo <- disa[ina == m, ]
              apo <- Rfast::colSort(apo)
              ta[, m] <- Rfast::colmeans( apo[1:k[j], , drop = FALSE] )
            }
            g[, j] <- Rfast::rowMins(ta)
          }
	      } else {
	        for (j in 1:klen) {
            for (m in 1:nc) {
              apo <- disa[ina == m, ]
              apo <- Rfast::colSort(apo)
              ta[, m] <- Rfast::colhameans( apo[1:k[j], , drop = FALSE] )
            }
		        g[, j] <- Rfast::rowMins(ta)
          }
	      }  ## end if (mesos)
      } else {
        for (i in 1:nu) {
          znewi <- znew[i, ]
          sa <- ( tzx - znewi )^2 / ( zx + znewi )
          sa[is.infinite(sa)] <- 0
          di <- Rfast::colsums( sa )
          di <- Rfast::Order(di)
          di <- di[ di <= max(k) ]
          for (j in 1:klen) {
            ind <- di[ 1:k[j] ]
            a <- Rfast::Table( ina[ind] )
            b <- as.numeric( names(a) )
            g[i, j] <- b[which.max(a)]
          }  ## end for (j in 1:klen)
        }  ## end for (i in 1:nu)
      }  ## end if (type == "NS")
      ## disa <- sqrt(disa) / abs(a) * sqrt(2 * p) not necessary to take the sqrt and then divide and multiply with constants everywhere
    }  ## end if (apostasi == "CS")
  }  ## end of other methods
  colnames(g) <- paste("k=", k, sep = "")
  g
}
