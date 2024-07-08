tflr.indeptest <- function(y, x, R = 999, ncores = 1) {

  kl <- Compositional::tflr(y, x)$kl
  n <- dim(y)[1]
  pkl <- numeric(R)

  if ( ncores <= 1 ) {
    for ( i in 1:R ) {
      id <- Rfast2::Sample.int(n, n)
      pkl[i] <- Compositional::tflr(y, x[id, ])$kl
    }

  } else {
    requireNamespace("doParallel", quietly = TRUE, warn.conflicts = FALSE)
    cl <- parallel::makePSOCKcluster(ncores)
    doParallel::registerDoParallel(cl)
    pkl <- foreach::foreach(i = 1:R, .combine = "c",
                            .packages = c("Compositional", "Rfast", "Rfast2") ) %dopar% {
                              id <- Rfast2::Sample.int(n, n)
                              return( Compositional::tflr(y, x[id, ])$kl )
                            }
  }

  ( sum(pkl < kl) + 1 ) / (R + 1)
}












# tflr.indeptest <- function(y, x, R = 999, ncores = 1) {
#   py <- dim(y)[2]   ;    px <- dim(x)[2]
#   pyx <- py * px    ;    n <- dim(y)[1]
#   n <- dim(y)[1]    ;    npy <- n * py
#
#   X <- matrix(0, npy, pyx)
#   indr <- matrix( 1:npy, ncol = py )
#   indc <- matrix( 1:pyx, ncol = py )
#   for ( i in 1:py )  X[ indr[, i], indc[, i] ] <- x
#   Y <- as.vector(y)
#
#   ind <- matrix( 1:pyx, ncol = px, byrow = TRUE )
#   A <- matrix(0, pyx, pyx)
#   for ( i in 1:px )  A[i, ind[, i]] <- 1
#   A <- t( rbind( A, diag(pyx), -diag(pyx) ) )
#   A <- A[, -c( (px + 1): pyx) ]
#   bvec <- c( rep(1, px), rep(0, pyx), rep(-1, pyx) )
#   A <- t(A)
#
#   a <- goric::orglm(Y ~ X - 1, family = quasibinomial(link="identity"),
#                     data = data.frame(Y = Y, X = X), constr = A, rhs = bvec, nec = px)
#   be <- matrix( abs(coef(a)), ncol = py)
#   est <- x %*% be
#   kl <- sum(y * log(y / est), na.rm = TRUE)
#
#   pkl <- rep(NA, R)
#   if ( ncores <= 1 ) {
#     for (i in 1:R) {
#       id <- Rfast2::Sample.int(n,n)
#       yp <- y[id, ]
#       Yp <- as.vector(yp)
#       a <- goric::orglm(Yp ~ X - 1, family = quasibinomial(link="identity"),
#                         data = data.frame(Yp = Yp, X = X), constr = A, rhs = bvec, nec = px)
#       be <- matrix( abs(coef(a)), ncol = py)
#       est <- x %*% be
#       pkl[i] <- sum(yp * log(yp / est), na.rm = TRUE)
#     }
#
#   } else {
#     requireNamespace("doParallel", quietly = TRUE, warn.conflicts = FALSE)
#     cl <- parallel::makePSOCKcluster(ncores)
#     doParallel::registerDoParallel(cl)
#     pkl <- foreach::foreach(i = 1:R, .combine = "c",
#                             .packages = c("Compositional", "Rfast", "Rfast2", "goric") ) %dopar% {
#                               id <- Rfast2::Sample.int(n, n)
#                               yp <- y[id, ]
#                               Yp <- as.vector(yp)
#                               a <- goric::orglm(Yp ~ X - 1, family = quasibinomial(link="identity"),
#                                                 data = data.frame(Yp = Yp, X = X), constr = A, rhs = bvec, nec = px)
#                               be <- matrix( abs(coef(a)), ncol = py)
#                               est <- x %*% be
#                               return( sum(yp * log(yp / est), na.rm = TRUE) )
#                             }
#   }
#
#   ( sum(pkl < kl) + 1 ) / (R + 1)
# }








