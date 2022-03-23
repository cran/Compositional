################################
#### Estimating parameters for compositional data
#### using the additive or the isometric log-ratio transformation
#### Tsagris Michail 1/2016
#### mtsagris@yahoo.gr
################################
 comp.den <- function(x, type = "alr", dist = "normal", tol = 1e-07) {
  if (dist == "normal") {
    if (type == "alr") {  ## additive log-ratio transformation
      y <- Compositional::alr(x)
      m <- Rfast::colmeans(y)
      mu <- c( 1, exp(m) )
      mu <- mu / sum(mu)
      s <- Rfast::cova(y)
    } else {
      y <- alfa(x, 0)$aff
      m <- Rfast::colmeans(y)
      mu <- Compositional::alfainv(m, 0)
      s <- Rfast::cova(y)
    }
    result <- list(mean = m, comp.mean = mu, covariance = s)

  } else if (dist == "t") {
      if (type == "alr") {  ## additive log-ratio transformation
        y <- Compositional::alr(x)
        mod <- Compositional::multivt(y)
        m <- mod$center
        mu <- c( 1, exp(m) )
        mu <- mu / sum(mu)
        s <- mod$scatter
        dof <- mod$dof
      } else {
        y <- Compositional::alfa(x, 0)$aff
        mod <- Compositional::multivt(y)
        m <- mod$center
        mu <- alfainv(m, 0)
        s <- mod$scatter
        dof <- mod$dof
      }
      result <- list(mean = m, comp.mean = mu, covariance = s, dof = dof)

  } else if (dist == "rob") {
      if (type == "alr") {  ## additive log-ratio transformation
        y <- Compositional::alr(x)
        mod <- MASS::cov.rob(y, method = "mcd")
        m <- mod$center
        mu <- c( 1, exp(m) )
        mu <- mu / sum(mu)
        s <- mod$cov
        best <- mod$best
      } else {
        y <- Compositional::alfa(x, 0)$aff
        mod <- MASS::cov.rob(y, method = "mcd")
        m <- mod$center
        mu <- alfainv(m, 0)
        s <- mod$cov
        best <- mod$best
      }
      result <- list(mean = m, comp.mean = mu, covariance = s, best = best)

  } else if (dist == "spatial") {
      if (type == "alr") {  ## additive log-ratio transformation
        y <- Compositional::alr(x)
        delta <- Rfast::spat.med(y, tol = tol)
        comp.delta <- c( 1, exp( delta ) )
        comp.delta <- delta / sum( delta )
        s <- Rfast::sscov(y, delta)
      } else {
        y <- Compositional::alfa(x, 0)$aff
        delta <- Rfast::spat.med(y)
        comp.delta <- alfainv(delta, 0)
        s <- Rfast::sscov(y, delta)
      }
      result <- list(spatmed = delta, comp.spat.med = comp.delta, ssc = s)

  } else if (dist == "skewnorm") {
      if (type == "alr") {  ## additive log-ratio transformation
        y <- Compositional::alr(x)
        mod <- sn::msn.mle(y = y)
        beta <- as.vector( mod$dp$beta )
        Omega <- as.matrix( mod$dp$Omega )
        alpha <- as.vector(mod$dp$alpha)
        cobeta <- c( 1, exp( beta) )
        cobeta <- cobeta / sum(cobeta)
      } else {
        y <- Compositional::alfa(x, 0)$aff
        mod <- sn::msn.mle(y = y)
        beta <- as.vector( mod$dp$beta )
        Omega <- as.matrix( mod$dp$Omega )
        alpha <- as.vector(mod$dp$alpha)
        cobeta <- Compositional::alfainv(beta, 0)
      }
      result <- list(beta = beta, Omega = Omega, alpha = alpha, comp.beta = cobeta)
  }

  result
}
