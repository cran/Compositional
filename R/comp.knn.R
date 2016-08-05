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
  x <- x / as.vector( Rfast::rowsums(x) )  ## makes sure the data sum to 1
  n <- nrow(x)
  p <- ncol(x)
  ina <- as.numeric(ina)
  xnew <- as.matrix(xnew)
  xnew <- matrix( xnew, ncol = p ) ## makes sure xnew is a matrix
  xnew <- xnew / as.vector( Rfast::rowsums(xnew) )  ## make the data sum to 1
  nc <- max(ina)  ## The number of groups
  nu <- nrow(xnew)
  disa <- matrix(0, nu, n)

  if (apostasi == "CS" & a == 0) {
    apostasi = "Ait"
  }

  if (apostasi == "ESOV") {
    xa <- x^a
    zx <- xa / as.vector( Rfast::rowsums( xa ) ) ## The power transformation is applied
    za <- xnew^a
    znew <- za / as.vector( Rfast::rowsums( za ) ) ## The power transformation is applied

    for (i in 1:nu) {
      zan <- znew[i, ]
      for (j in 1:n) {
        ma <- zan + zx[j, ]
        disa[i, j] <- sqrt( sum( zan * log( 2 * zan/ma ) +
                      zx[j, ] * log( 2 * zx[j, ]/ma ), na.rm = TRUE ) )
      }
    }

  } else  if ( apostasi == "taxicab" ) {
    xa <- x^a
    zx <- xa / as.vector( Rfast::rowsums( xa ) ) ## The power transformation is applied
    za <- xnew^a
    znew <- za / as.vector( Rfast::rowsums( za ) ) ## The power transformation is applied

    for (i in 1:nu) {
      b <- t(zx) - znew[i, ]
      disa[i, ] <- as.vector( Rfast::colsums( abs(b) ) )
    }

  } else if ( apostasi == "Ait" ) {
    ## this requires non zero data ## be careful
    xa <- log(x)
    zx <- xa - as.vector( Rfast::rowmeans( xa ) )
    za <- log(xnew)
    znew <- za - as.vector( Rfast::rowmeans( za ) )
    tzx <- t(zx)

    for (i in 1:nu) {
      zz <- tzx - znew[i, ]
      disa[i, ] <-  sqrt( as.vector( Rfast::colsums( zz^2 ) ) )
    }

  } else if ( apostasi == "Hellinger" ) {
    zx <- sqrt(x)
    znew <- sqrt(xnew)
    tzx <- t(zx)

    for (i in 1:nu) {
      zz <- tzx - znew[i, ]
      disa[i, ] <-  sqrt( as.vector( Rfast::colsums( zz^2 ) ) )
    }
    disa <- disa / sqrt(2)

  } else if ( apostasi == "angular" ) {
    zx <- sqrt(x)
    znew <- sqrt(xnew)
    disa <- tcrossprod( znew, zx )
    disa[disa >= 1] <- 1
    disa <- acos(disa)

  }  else if ( apostasi == "CS" ) {
     xa <- x^a
     zx <- xa / as.vector( Rfast::rowsums( xa ) ) ## The power transformation is applied
     za <- xnew^a
    znew <- za / as.vector( Rfast::rowsums( za ) ) ## The power transformation is applied

     for (i in 1:nu) {
       for (j in 1:n) {
         sa <- ( zx[j, ] - znew[i, ] )^2 / (zx[j, ] + znew[i, ])
         disa[i, j] <- sum( sa[ abs(sa)<Inf ] )
       }
     }
     disa <- 1 / abs(a) * sqrt(2 * p) * sqrt(disa)
  }

  ta <- matrix(nrow = nu, ncol = nc)

  if (type == "NS") {
    ## Non Standard algorithm

    for (m in 1:nc) {
      disa <- t( apply(disa, 1, sort) )
      if (mesos == TRUE) {
        ta[, m] <- as.vector( Rfast::rowmeans( disa[, 1:k] ) )

      } else {
        ta[, m] <- k / as.vector( Rfast::rowsums( 1 / disa[, 1:k] ) )
      }
    }

    g <- apply(ta, 1, which.min)

  } else {   ## if type is "S"
    ## Standard algorithm

    g <- numeric(nu)
    for (l in 1:nu) {
      xa <- cbind(ina, disa[l, ])
      qan <- xa[order(xa[, 2]), ]
      sa <- qan[1:k, 1]
      tab <- table(sa)
      g[l] <- as.integer(names(tab)[ which.max(tab) ] )
    }
  }

  g
}
