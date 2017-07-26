comp.knn <- function(xnew, x, ina, a = 1, k = 5, type = "S", apostasi = "ESOV", mesos = TRUE) {
  n <- dim(x)[1]
  p <- dim(x)[2]
  ina <- as.numeric(ina)
  xnew <- as.matrix(xnew)
  xnew <- matrix( xnew, ncol = p ) ## makes sure xnew is a matrix
  nc <- max(ina)  ## The number of groups
  nu <- nrow(xnew)

  if (apostasi == "CS" & a == 0)  apostasi = "Ait"

  if ( (apostasi == "taxicab" | apostasi == "Ait" | apostasi == "Hellinger")  &  type == "S" ) {
    if ( mesos ) {
      method <- "average"
    } else  method <- "median"

    if ( apostasi == "taxicab" ) {
      xa <- x^a
      zx <- xa / Rfast::rowsums( xa )  ## The power transformation is applied
      za <- xnew^a
      znew <- za / Rfast::rowsums( za )  ## The power transformation is applied
      g <- Rfast::knn(znew, ina, zx, k = k, dist.type = "mahattan", type = "C",
                      method = method, freq.option = 1)
    } else if ( apostasi == "Ait" ) {
      xa <- Rfast::Log(x)
      zx <- xa - Rfast::rowmeans( xa )
      za <- Rfast::Log(xnew)
      znew <- za - Rfast::rowmeans( za )
      g <- Rfast::knn(znew, ina, zx, k = k, dist.type = "euclidean", type = "C",
                      method = method, freq.option = 1)

    } else if ( apostasi == "Hellinger" ) {
      g <- Rfast::knn(sqrt(xnew), ina, sqrt(x), k = k, dist.type = "euclidean", type = "C",
                      method = method, freq.option = 1)
    }

  } else {

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
      disa <- Rfast::dista(sqrt(xnew), sqrt(x), "manhattan", trans = FALSE)
    }

    if (apostasi == "ESOV") {
      disa <- matrix(0, n, nu)
      xa <- x^a
      zx <- xa / Rfast::rowsums( xa )  ## The power transformation is applied
      za <- xnew^a
      znew <- za / Rfast::rowsums( za )  ## The power transformation is applied
      for (i in 1:nu) {
        zan <- znew[i, ]
        for (j in 1:n) {
          zxj <- zx[j, ]
          ma <- zan + zxj
          disa[j, i] <- sqrt( sum( zan * log( 2 * zan / ma ) + zxj * log( 2 * zxj/ma ), na.rm = TRUE ) )
        }
      }
    } else if ( apostasi == "angular" ) {
      zx <- sqrt(x)
      znew <- sqrt(xnew)
      disa <- tcrossprod(zx, znew )
      disa[disa >= 1] <- 1
      disa <- acos(disa)

    } else if ( apostasi == "CS" ) {
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
      disa <- sqrt(disa) ## / abs(a) * sqrt(2 * p) not necessary to divide and multiply with constants everywhere
    }

    if (type == "NS") {      ## Non Standard algorithm
      ta <- matrix(nrow = nu, ncol = nc)

      for (m in 1:nc) {
        apo <- disa[ina == m, ]
        apo <- Rfast::sort_mat(apo)
        if ( mesos ) {
          ta[, m] <- Rfast::colmeans( apo[1:k, ] )
        } else  ta[, m] <- k / Rfast::colsums( 1 / apo[1:k, ] )
      }
      g <- Rfast::rowMins(ta)

    } else {   ## if type is "S"   ## Standard algorithm
      g <- numeric(nu)
      for (l in 1:nu) {
        xa <- cbind(ina, disa[, l])
        qan <- xa[order(xa[, 2]), ]
        sa <- qan[1:k, 1]
        tab <- table(sa)
        g[l] <- as.integer(names(tab)[ which.max(tab) ] )
      }
    }

  }
  g
}
