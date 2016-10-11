################################
#### Classification for compositional data using a power transformation
#### The k-NN algorithm
#### Tsagris Michail 7/2015
#### References: Tsagris, M. T. (2014).
#### The k-NN algorithm for compositional data: a revised approach with and without zero values present
#### Journal of Data Science, 12(3):519-534
#### mtsagris@yahoo.gr
################################

comp.knn <- function(xnew, x, ina, a = 1, k = 5, type = "S",
                     apostasi = "ESOV", mesos = TRUE) {
  ## xnew is the new dataset. It can be a single vector or a matrix
  ## x is the matrix containing the data
  ## ina indicates the groups
  ## a is the value of the power parameter
  ## k in the number of nearest neighbours
  ## apostasi is the type of metric used, "ESOV", "taxicab",
  ## "Ait", "Hellinger", "angular" or "CS"
  ## type is either S or NS. Should the standard k-NN be use or not
  ## if mesos is TRUE, then the arithmetic mean distance of the
  ## k nearest points will be used
  ## If not, then the harmonic mean will be used.
  ## Both of these apply for the non-standard,
  ## algorithm, that is when type=NS

  x <- as.matrix(x)  ## makes sure x is a matrix
  x <- x / Rfast::rowsums(x)  ## makes sure the data sum to 1
  n <- dim(x)[1]
  p <- dim(x)[2]
  ina <- as.numeric(ina)
  xnew <- as.matrix(xnew)
  xnew <- matrix( xnew, ncol = p ) ## makes sure xnew is a matrix
  xnew <- xnew / Rfast::rowsums(xnew)  ## make the data sum to 1
  nc <- max(ina)  ## The number of groups
  nu <- nrow(xnew)
  disa <- matrix(0, n, nu)

  if (apostasi == "CS" & a == 0) {
    apostasi = "Ait"
  }

  if (apostasi == "ESOV") {
    xa <- x^a
    zx <- xa / Rfast::rowsums( xa )  ## The power transformation is applied
    za <- xnew^a
    znew <- za / Rfast::rowsums( za )  ## The power transformation is applied

    for (i in 1:nu) {
      zan <- znew[i, ]
      for (j in 1:n) {
        zxj <- zx[j, ]
        ma <- zan + zxj
        disa[j, i] <- sqrt( sum( zan * log( 2 * zan / ma ) +
                                   zxj * log( 2 * zxj/ma ), na.rm = TRUE ) )
      }
    }

  } else  if ( apostasi == "taxicab" ) {
    xa <- x^a
    zx <- xa / Rfast::rowsums( xa )  ## The power transformation is applied
    za <- xnew^a
    znew <- za / Rfast::rowsums( za )  ## The power transformation is applied

    for (i in 1:nu) {
      b <- t(zx) - znew[i, ]
      disa[, i] <- Rfast::colsums( abs(b) )
    }

  } else if ( apostasi == "Ait" ) {
    ## this requires non zero data ## be careful
    xa <- log(x)
    zx <- xa - Rfast::rowmeans( xa )
    za <- log(xnew)
    znew <- za - Rfast::rowmeans( za )
    tzx <- t(zx)

    for (i in 1:nu) {
      zz <- tzx - znew[i, ]
      disa[, i] <-  sqrt( Rfast::colsums( zz^2 ) )
    }

  } else if ( apostasi == "Hellinger" ) {
    zx <- sqrt(x)
    znew <- sqrt(xnew)
    tzx <- t(zx)

    for (i in 1:nu) {
      zz <- tzx - znew[i, ]
      disa[, i] <-  sqrt( Rfast::colsums( zz^2 ) )
    }
    disa <- disa / sqrt(2)

  } else if ( apostasi == "angular" ) {
    zx <- sqrt(x)
    znew <- sqrt(xnew)
    disa <- tcrossprod(zx, znew )
    disa[disa >= 1] <- 1
    disa <- acos(disa)

  }  else if ( apostasi == "CS" ) {
    xa <- x^a
    zx <- xa / Rfast::rowsums( xa )  ## The power transformation is applied
    za <- xnew^a
    znew <- za / Rfast::rowsums( za )  ## The power transformation is applied

    for (i in 1:nu) {
      znewi <- znew[i, ]
      for (j in 1:n) {
        zxj <- zx[j, ]
        sa <- ( zxj - znewi )^2 / ( zxj + znewi )
        disa[j, i] <- sum( sa[ abs(sa)<Inf ] )
      }
    }
    disa <- sqrt(2 * p) * sqrt(disa) / abs(a)
  }

  if (type == "NS") {
    ## Non Standard algorithm

    ta <- matrix(nrow = nu, ncol = nc)

    for (m in 1:nc) {
      apo <- disa[ina == m, ]
      apo <- Rfast::sort_mat(apo)
      if (mesos == TRUE) {
        ta[, m] <- Rfast::colmeans( apo[1:k, ] )
      } else {
        ta[, m] <- k / Rfast::colsums( 1 / apo[1:k, ] )
      }
    }

    g <- max.col(-ta)

  } else {   ## if type is "S"
    ## Standard algorithm

    g <- numeric(nu)
    for (l in 1:nu) {
      xa <- cbind(ina, disa[, l])
      qan <- xa[order(xa[, 2]), ]
      sa <- qan[1:k, 1]
      tab <- table(sa)
      g[l] <- as.integer(names(tab)[ which.max(tab) ] )
    }
  }

  g
}
