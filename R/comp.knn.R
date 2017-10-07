comp.knn <- function(xnew, x, ina, a = 1, k = 5, type = "S", apostasi = "ESOV", mesos = TRUE) {
  n <- dim(x)[1]
  p <- dim(x)[2]
  ina <- as.numeric(ina)
  xnew <- as.matrix(xnew)
  xnew <- matrix(xnew, ncol = p ) ## makes sure xnew is a matrix
  nc <- max(ina)  ## The number of groups
  nu <- dim(xnew)[1]

  if (apostasi == "CS" & a == 0)  apostasi = "Ait"

  if ( (apostasi == "taxicab" | apostasi == "Ait" | apostasi == "Hellinger")  &  type == "S" ) {

    if ( apostasi == "taxicab" ) {
      xa <- x^a
      zx <- xa / Rfast::rowsums( xa )  ## The power transformation is applied
      za <- xnew^a
      znew <- za / Rfast::rowsums( za )  ## The power transformation is applied
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
    if ( apostasi == "taxicab" ) {
      xa <- x^a
      zx <- xa / Rfast::rowsums( xa )  ## The power transformation is applied
      za <- xnew^a
      znew <- za / Rfast::rowsums( za )  ## The power transformation is applied
      disa <- Rfast::dista(znew, zx, "manhattan", trans = FALSE)

    } else if ( apostasi == "Ait" ) {
      xa <- Rfast::Log(x)
      zx <- xa - Rfast::rowmeans( xa )
      za <- Rfast::Log(xnew)
      znew <- za - Rfast::rowmeans( za )
      disa <- Rfast::dista(znew, zx, trans = FALSE)

    } else if ( apostasi == "Hellinger" ) {
      disa <- Rfast::dista(sqrt(xnew), sqrt(x), "euclidean", trans = FALSE)

    } else if ( apostasi == "angular" ) {
      zx <- sqrt(x)
      znew <- sqrt(xnew)
      disa <- tcrossprod(zx, znew )
      disa[disa >= 1] <- 1
      disa <- acos(disa)

    } else if ( apostasi == "ESOV" ) {
      disa <- matrix(0, n, nu)
      xa <- x^a
      zx <- xa / Rfast::rowsums( xa )  ## The power transformation is applied
      za <- xnew^a
      znew <- za / Rfast::rowsums( za )  ## The power transformation is applied
      tzx <- t(zx)
      for (i in 1:nu) {
        zan <- znew[i, ]
        ma <- tzx + zan
        disa[, i] <- colSums( zan * log( 2 * zan / ma ) + tzx * log( 2 * tzx/ma ), na.rm = TRUE )
      }
    } else if ( apostasi == "CS" ) {
      xa <- x^a
      zx <- xa / Rfast::rowsums( xa )  ## The power transformation is applied
      za <- xnew^a
      znew <- za / Rfast::rowsums( za )  ## The power transformation is applied
      tzx <- t(zx)
      for (i in 1:nu) {
        znewi <- znew[i, ]
        sa <- ( tzx - znewi )^2 / ( zx + znewi )
        sa[is.infinite(sa)] <- 0
        disa[, i] <- Rfast::colsums( sa )
      }
      ## disa <- sqrt(disa) / abs(a) * sqrt(2 * p) not necessary to take the sqrt and then divide and multiply with constants everywhere
    }

    klen <- length(k)
    if ( klen == 1 ) {

      if (type == "NS") {      ## Non Standard algorithm
        ta <- matrix(nrow = nu, ncol = nc)
        for (m in 1:nc) {
          apo <- disa[ina == m, ]
          apo <- Rfast::sort_mat(apo)
          if ( mesos ) {
            ta[, m] <- Rfast::colmeans( apo[1:k, , drop = FALSE] )
          } else  ta[, m] <- Rfast::colhameans( apo[1:k, , drop = FALSE] )
        }
        g <- as.matrix( Rfast::rowMins(ta) )

      } else {   ## if type is "S"   ## Standard algorithm
        g1 <- Rfast::colnth( disa, rep(k, nu) )
        g <- g1
        for (l in 1:nu) {
          ind <- which(disa[, l] <= g1[l] )
          a <- Rfast::Table( ina[ind] )
          b <- as.numeric( names(a) )
          g[l] <- b[which.max(a)]
        }
        g <- as.matrix(g)
      }  ## end if (type == "NS")

    } else {  ## k has many values
      g <- matrix(0, nu, klen)
      if (type == "NS") {      ## Non Standard algorithm
        ta <- matrix(nrow = nu, ncol = nc)
        for (j in 1:klen) {
          for (m in 1:nc) {
            apo <- disa[ina == m, ]
            apo <- Rfast::sort_mat(apo)
            if ( mesos ) {
              ta[, m] <- Rfast::colmeans( apo[1:k[j], , drop = FALSE] )
            } else  ta[, m] <- Rfast::colhameans( apo[1:k[j], , drop = FALSE] )
          }
          g[, j] <- Rfast::rowMins(ta)
        }

      } else {   ## if type is "S"   ## Standard algorithm
        for (j in 1:klen) {
          g1 <- Rfast::colnth( disa, rep( k[j], nu) )
          for (l in 1:nu) {
            ind <- which(disa[, l] <= g1[l] )
            a <- Rfast::Table( ina[ind] )
            b <- as.numeric( names(a) )
            g[l, j] <- b[which.max(a)]
          }  ## end inner for
        } ## end outer for
      }  ## end if (type == "NS")
    }  ## end if (length(k) == 1)
  }  ## end of other methods
  colnames(g) <- paste("k=", k, sep = "")
  g
}
