################################
#### Classification for compositional data using the alpha-transformation
#### The k-NN algorithm
#### Tsagris Michail 8/2015
#### References: Tsagris, M., Preston S. and Wood A.T.A. (2016).
#### Improved classication for compositional data using the alpha-transformation
#### Journal of Classification (To appear)
#### http://arxiv.org/pdf/1506.04976v2.pdf
#### mtsagris@yahoo.gr
################################

alfa.knn <- function(xnew, x, ina, a = 1, k = 5, type = "S", mesos = TRUE) {
  ## x is the matrix containing the data
  ## ina indicates the groups
  ## a is the value of the power parameter
  ## k in the number of nearest neighbours
  ## apostasi is the type of metric used, either ESOV_a or Taxicab_a
  ## type is either S or NS. Should the standard k-NN be use or not
  ## if mesos is TRUE, then the arithmetic mean distance of the
  ## k nearest points will be used
  ## If not, then the harmonic mean will be used.
  ## Both of these apply for the non-standard,
  ## algorithm, that is when type=NS
  ## xnew is the new dataset. It can be a single vector or a matrix

  x <- as.matrix(x)  ## makes sure x is a matrix
  x <- x / Rfast::rowsums(x)  ## makes sure the data sum to 1
  n <- dim(x)[1]
  p <- dim(x)[2]
  xnew <- as.matrix(xnew)
  xnew <- matrix(xnew, ncol = p)  ## makes sure xnew is a matrix
  xnew <- xnew / Rfast::rowsums(xnew)  ## make the data sum to 1
  ina <- as.numeric(ina)
  nc <- max(ina) ## The number of groups
  nu <- nrow(xnew)
  apo <- matrix( 0, n, nu )

  znew <- alfa(xnew, a)$aff
  z <- alfa(x, a)$aff
  tz <- t(z)

  for (i in 1:nu) {
    zz <- tz - znew[i, ]
    apo[, i] <- sqrt( Rfast::colsums( zz^2 ) )
  }

  if (type == "NS") {
    ## Non Standard algorithm
    ta <- matrix(nrow = nu, ncol = nc)
    for (m in 1:nc) {
      dista <- apo[ina == m, ]
      dista <- Rfast::sort_mat(dista)
      if (mesos == TRUE) {
        ta[, m] <- Rfast::colmeans( dista[1:k, ] )
      } else {
        ta[, m] <- k / Rfast::colsums( 1 / dista[1:k, ] )
      }
    }
    g <- max.col(-ta)

  } else if (type == "S") {
    ## Standard algorithm
    g <- numeric(nu)
    for (l in 1:nu) {
      xa <- cbind(ina, apo[, l])
      qan <- xa[order(xa[, 2]), ]
      sa <- qan[1:k, 1]
      tab <- table(sa)
      g[l] <- as.integer( names(tab)[ which.max(tab) ] )
    }
  }

  g
}
