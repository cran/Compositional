rda <- function(xnew, x, ina, gam = 1, del = 0) {
  n <- dim(x)[1]
  D <- dim(x)[2]
  xnew <- as.matrix(xnew)
  xnew <- matrix(xnew, ncol = D)
  nu <- dim(xnew)[1]  ## number of the new observations
  ina <- as.numeric(ina)
  ng <- tabulate(ina)
  nc <- length(ng)
  ta <- matrix(nrow = nu, ncol = nc)
  ci <- 2 * log(ng / n)
  sk <- vector("list", nc)
  mesos <- Rfast2::colGroup(x, as.integer(ina) ) / ng
  sa <- 0
  
  for (i in 1:nc) {
    xi <- x[ina == i, ]
    m <- sqrt(ng[i]) * mesos[i, ]
    sk[[ i ]] <- ( crossprod(xi) - tcrossprod(m) )
    sa <- sa + sk[[ i ]]
    sk[[ i ]] <- sk[[ i ]] / (ng[i] - 1)   
  }
  
  Sp <- sa/(n - nc)
  sp <- diag( sum( diag( Sp ) ) / D, D ) ## spherical covariance matrix
  Sa <- gam * Sp + (1 - gam) * sp  ## regularised covariance matrix

  for (j in 1:nc) {
    Ska <- del * sk[[ j ]] + (1 - del) * Sa
    ta[, j] <- ci[j] - log( det( Ska ) ) - Rfast::mahala( xnew, mesos[j, ], Ska )
    ## the scores are doubled, but for efficiency I did not multiply with 0.5
  }

  est <- Rfast::rowMaxs(ta)
  expta <- exp(ta)
  prob <- expta / Rfast::rowsums( expta ) ## the probability of classification
  list(prob = prob, scores = ta, est = est)
}
