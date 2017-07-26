compknn.tune <- function(x, ina, M = 10, A = 5, type = "S", mesos = TRUE, a = seq(-1, 1, by = 0.1), apostasi = "ESOV", mat = NULL, graph = FALSE) {
  n <- dim(x)[1]  ## sample size
  ina <- as.numeric(ina)
  if ( A >= min(table(ina)) )  A <- min( table(ina) ) - 3  ## The maximum
  ## number  of nearest neighbours to use
  ng <- max(ina)  ## The number of groups
  if ( min(x) == 0 )  a <- a[ a > 0 ]
  dis <- matrix(0, n, n)
  if ( is.null(mat) ) {
    nu <- sample(1:n, min( n, round(n / M) * M ) )
    ## It may be the case this new nu is not exactly the same
    ## as the one specified by the user
    ## to a matrix a warning message should appear
    options(warn = -1)
    mat <- matrix( nu, ncol = M ) # if the length of nu does not fit
  } else  mat <- mat

  M <- dim(mat)[2]
  rmat <- dim(mat)[1]
  ## The algorithm is repated R times and each time the estimated
  ## percentages are stored in the array per.
  if (apostasi == "ESOV" | apostasi == "taxicab" | apostasi == "CS") {

    runtime <- proc.time()
    a <- a[ a != 0 ]

    per <- array( dim = c(M, A - 1, length(a)) )

    for ( i in 1:length(a) ) {

      z <- x^a[i] / Rfast::rowsums( x^a[i] )  ## The power transformation is applied

      if (apostasi == "ESOV") {
        for ( m1 in 1:c(n - 1) ) {
          z1 <- z[m1, ]
          for ( m2 in c(m1 + 1):n ) {
            z2 <- z[m2, ]
            ma <- z1 + z2
            dis[m1, m2] <- sqrt( sum( z1 * log( 2 * z1 / ma ) + z2 * log( 2 * z2 / ma ), na.rm = TRUE ) )
          }
        }
        dis <- dis + t(dis)

      } else if (apostasi == "taxicab") {
        dis <- Rfast::Dist(z, method = "manhattan")

      } else if ( apostasi == "CS" ) {
        p <- dim(x)[2]
        for ( m1 in 1:c(n - 1) ) {
          z1 <- z[m1, ]
          for ( m2 in c(m1 + 1):n ) {
            z2 <- z[m2, ]
            sa <- (z1 - z2)^2 / (z1 + z2)
            dis[m1, m2] <- sum( sa[ abs(sa) < Inf ] )
          }
        }
        dis <- sqrt(dis) / abs( a[i] )  ## * sqrt(2 * p)  not necessary to multiply with constant everything
        dis <- dis + t(dis)
      }
      ## The k-NN algorithm is calculated R times. For every repetition a
      ## test sample is chosen and its observations are classified
      for (vim in 1:M) {
        id <- ina[ mat[, vim] ]  ## groups of test sample
        ina2 <- ina[ -mat[, vim] ]   ## groups of training sample
        aba <- as.vector( mat[, vim] )
        aba <- aba[aba > 0]
        apo <- dis[-aba, aba]

        if (type == "NS") {
          ## Non Standard algorithm
          ta <- matrix(nrow = rmat, ncol = ng)
          for ( j in 1:c(A - 1) ) {
            knn <- j + 1
            for (l in 1:ng) {
              dista <- apo[ina2 == l, ]
              dista <- Rfast::sort_mat(dista)
              if ( mesos ) {
                ta[, l] <- Rfast::colmeans( dista[1:knn, ] )
              } else  ta[, l] <- knn / Rfast::colsums( 1 / dista[1:knn, ] )
            }
            g <- Rfast::rowMins(ta)
            per[vim, j, i] <- sum( g == id ) / rmat
          }

        } else if (type == "S") {
          ## Standard algorithm
          for (j in 1:c(A - 1) ) {
            g <- numeric(rmat)
            knn <- j + 1
            for (k in 1:rmat) {
              xa <- cbind(ina2, apo[, k])
              qan <- xa[order(xa[, 2]), ]
              sa <- qan[1:knn, 1]
              tab <- table(sa)
              g[k] <- as.integer(names(tab)[which.max(tab)])
            }
            per[vim, j, i] <- sum( g == id ) / rmat
          }
        }
      }
    }

    ela <- matrix(nrow = length(a), ncol = A - 1)
    for ( i in 1:length(a) )  ela[i, ] <- colMeans(per[, , i])
    ## The ela matrix contains the averages of the R
    ## repetitions over alpha and k
    colnames(ela) <- paste("k=", 2:A, sep = "")
    rownames(ela) <- paste("alpha=", a, sep = "")
    ## The code for the heat plot of the estimated percentages
    if (graph)  fields::image.plot(a, 2:A, ela, col = grey(1:11/11), ylab = "k nearest-neighbours", xlab = expression(paste(alpha, " values")) )

    opt <- max(ela)
    confa <- which(ela == opt, arr.ind = TRUE)[1, ]
    bias <- numeric(M)
    for (i in 1:M)  bias[i] <- opt - per[ i, confa[2], confa[1] ]
    bias <- mean(bias)
    performance <- c(opt - bias, bias)
    names(performance) <- c( "rate", "bias" )
    runtime <- proc.time() - runtime
    results <- list( ela = ela, performance = performance, best_a = a[ confa[1] ], best_k = confa[2] + 1, runtime = runtime )

  } else if (apostasi == "Ait" | apostasi == "Hellinger" | apostasi == "angular" ) {

    runtime <- proc.time()
    per <- matrix(nrow = M, ncol = A - 1)

    if (apostasi == "Ait") {
      xa <- Rfast::Log(x)
      z <- xa - Rfast::rowmeans( xa )
      dis <- Rfast::Dist(z)
    } else if (apostasi == "Hellinger") {
      dis <- Rfast::Dist(z, "hellinger")
    } else if (apostasi == "angular") {
      dis <- tcrossprod( sqrt(x) )
      diag(dis) <- 1
      dis[ dis > 1 ] <- 1
      dis <- acos(dis)
    }
    diag(dis) <- 0

    for (vim in 1:M) {
      id <- ina[ mat[, vim] ]  ## groups of test sample
      ina2 <- ina[ -mat[, vim] ]   ## groups of training sample
      aba <- as.vector( mat[, vim] )
      aba <- aba[aba > 0]
      apo <- dis[-aba, aba]
      ta <- matrix(nrow = rmat, ncol = ng)

      if (type == "NS") {
        ## Non Standard algorithm
        for ( j in 1:c(A - 1) ) {
          knn <- j + 1
          for (l in 1:ng) {
            dista <- apo[ina2 == l, ]
            dista <- Rfast::sort_mat(dista)
            if ( mesos ) {
              ta[, l] <- Rfast::colmeans( dista[1:knn, ] )
            } else  ta[, l] <- knn / Rfast::colsums( 1 / dista[1:knn, ] )
          }
          g <- Rfast::rowMins(ta)
          per[vim, j] <- sum( g == id )/rmat
        }

      } else {   ## if (type == "S")
        ## Standard algorithm
        for ( j in 1:c(A - 1) ) {
          knn <- j + 1
          g <- numeric(rmat)
          for (k in 1:rmat) {
            xa <- cbind(ina2, apo[, k])
            qan <- xa[order(xa[, 2]), ]
            sa <- qan[1:knn, 1]
            tab <- table(sa)
            g[k] <- as.integer(names(tab)[which.max(tab)])
          }
          per[vim, j] <- sum( g == id )/rmat
        }
      }
    }

    ela <- Rfast::colmeans(per)
    opt <- max(ela)
    names(ela) <- paste("k=", 2:A, sep = "")
    best_k <- which.max(ela) + 1
    bias <- Rfast::rowMaxs(per, value = TRUE) - per[, best_k]  ## apply(per, 1, max) - per[, best_k]
    bias <- mean(bias)
    performance <- c(opt - bias, bias)
    names(performance) <- c( "rate", "bias" )

    if (graph)  plot(2:A, ela, type = "b", xlab = "k nearest neighbours", pch = 9, col = 2, ylab = "Estimated percentage of correct classification")
    runtime <- proc.time() - runtime
    results <- list(ela = ela, performance = performance, best_k = which.max(ela) + 1, runtime = runtime)
  }
  results
}
