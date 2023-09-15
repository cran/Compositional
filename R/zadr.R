zadr <- function(y, x, con = TRUE, B = 1, ncores = 2, xnew = NULL) {
  ## y is the compositional data
  ## x is the independent variable(s)
  dm <- dim(y)
  D <- dm[2]   ;   d <- D - 1
  ## d is the dimensionality of the simplex
  n <- dm[1] ## sample size

  beta.ini <- Compositional::kl.compreg(y, x, con = con)$be
  ini.phi <- Compositional::zad.est(y)$phi
  xini <- x

  x <- model.matrix(~., as.data.frame(x) )
  if ( !con )  x <- x[, -1, drop = FALSE]
  runtime <- proc.time()
  ## next we separate the compositional vectors, those which contain
  ## zeros and those without. The same separation is performed for the
  ## independent variable(s)
  a1 <- which( Rfast::rowsums( y > 0 ) == D )
  a2 <- which( Rfast::rowsums( y > 0 ) != D )
  n1 <- length(a1)
  n2 <- n - n1
  ## n1 is the sample size of the compositional vectors with no zeros
  ## n2 is the sample size of the compositional vectors with zeros
  za <- y[a2, , drop = FALSE]
  za[za == 0] <- 1
  za[ za < 1 ] <- 0
  theta <- table( apply(za, 1, paste, collapse = ",") )
  theta <- as.vector(theta)
  const <- n1 * log(n1/n) + sum( theta * log(theta/n) )

  y1 <- y[a1, , drop = FALSE]
  ly1 <- log( y1 )
  x1 <- x[a1, , drop = FALSE]
  ly2 <- log( y[a2, , drop = FALSE] )
  x2 <- x[a2, , drop = FALSE]
  n1 <- nrow(y1)    ;    n2 <- n - n1

  ##############
  ini.par <- c( log(ini.phi), as.vector(beta.ini) )  ## initial parameter values
  z <- list(ly1 = ly1, ly2 = ly2, x1 = x1, x2 = x2, a1 = a1, a2 = a2)
  suppressWarnings({
    qa <- optim( ini.par, .mixreg, z = z )
    qa <- optim( qa$par, .mixreg, z = z )
    qa <- optim( qa$par, .mixreg, z = z, hessian = TRUE, control = list(maxit = 1000) )
  })
  phi <- exp( qa$par[1] )  ## final phi value
  be <- matrix( qa$par[-1], ncol = d ) ## final beta values

  if ( B > 1 ) {
    runtime <- proc.time()
    requireNamespace("doParallel", quietly = TRUE, warn.conflicts = FALSE)
    cl <- parallel::makePSOCKcluster(ncores)
    doParallel::registerDoParallel(cl)
    betaboot <- foreach::foreach( i = 1:B, .combine = rbind, .export = c("Sample.int", ".mixreg", "diri.nr"),
                                 .packages = c("Rfast2", "Compositional") ) %dopar% {
      ida <- Rfast2::Sample.int(n, n, replace = TRUE)
      yb <- y[ida, ]
      xb <- x[ida, ]

      a1 <- which( Rfast::rowsums( yb > 0 ) == D )
      a2 <- which( Rfast::rowsums( yb > 0 ) != D )
      n1 <- length(a1)
      n2 <- n - n1
      ## n1 is the sample size of the compositional vectors with no zeros
      ## n2 is the sample size of the compositional vectors with zeros
      za <- yb[a2, , drop = FALSE]
      za[za == 0] <- 1
      za[ za < 1 ] <- 0
      theta <- table( apply(za, 1, paste, collapse = ",") )
      theta <- as.vector(theta)
      const <- n1 * log(n1/n) + sum( theta * log(theta/n) )

      y1 <- yb[a1, , drop = FALSE]
      ly1 <- log( y1 )
      x1 <- xb[a1, , drop = FALSE]
      ly2 <- log( y[a2, , drop = FALSE] )
      x2 <- xb[a2, , drop = FALSE]
      n1 <- nrow(y1)    ;    n2 <- n - n1
      ##############
      beta.ini <- try( Compositional::kl.compreg(yb, xini[ida, ], con = con)$be, silent = TRUE )
      if ( identical(class(beta.ini), "try-error") ) {
        beta.ini <- be
      }
      ini.phi <- Compositional::zad.est(yb)$phi
      ini.par <- c( log(ini.phi), as.vector(beta.ini) )  ## initial parameter values
      z <- list(ly1 = ly1, ly2 = ly2, x1 = x1, x2 = x2, a1 = a1, a2 = a2)
      suppressWarnings({
        qa <- optim( ini.par, .mixreg, z = z )
        qa <- optim( qa$par, .mixreg, z = z )
        qa <- optim( qa$par, .mixreg, z = z, control = list(maxit = 1000) )
      })
      return( qa$par )
    }  ##  end foreach
    parallel::stopCluster(cl)
    sigma <- cov(betaboot)
    seb <- sqrt( diag(sigma) )
    seb <- matrix(seb[-1], ncol = d)
    colnames(seb) <- colnames( y[, -1] )
    rownames(seb) <- colnames(x)
    runtime <- proc.time() - runtime

  } else {
    sigma <- try( solve( qa$hessian ), silent = TRUE )
    if ( !identical( class(sigma), "try-error") ) {
      seb <- sqrt( diag(sigma) )
      seb <- matrix(seb[-1], ncol = d)
      colnames(seb) <- colnames( y[, -1] )
      rownames(seb) <- colnames(x)
    } else {
      sigma <- NULL
      seb <- NULL
    }
  }  ##  end if (B > 1)

  colnames(be) <- colnames( y[, -1] )
  rownames(be) <- colnames(x)
  est <- NULL

  if ( !is.null(xnew) ) {
    xnew <- model.matrix(~., as.data.frame(xnew) )
    if ( !con )  xnew <- xnew[, -1, drop = FALSE]
    ma <- cbind(1, exp( xnew %*% be ) )
    est <- ma / Rfast::rowsums(ma)  ## fitted values
    colnames(est) <- colnames(y)
  }
  runtime <- proc.time() - runtime

  list(runtime = runtime, loglik = -qa$value + const, phi = phi, be = be, seb = seb, sigma = sigma, est = est )
}




.mixreg <- function(param, z) {
  ## separation of phi and the betas
  ## and the exponential to avoid negative values of phi
  phi <- exp( param[1] )   ;   para <- param[-1]
  ly1 <- z$ly1     ;     ly2 <- z$ly2
  x1 <- z$x1     ;     x2 <- z$x2
  a1 <- z$a1     ;     a2 <- z$a2
  dm <- dim(ly1)
  D <- dm[2]    ;    d <- D - 1
  n1 <- length(a1)    ;   n2 <- length(a2)
  ## n1 is the sample size of the compositional vectors with no zeros
  ## n2 is the sample size of the compositional vectors with zeros
  n <- n1 + n2  ## total sample size
  ## next we separate the compositional vectors, those which contain
  ## zeros and those without. The same separation is performed for the
  ## independent variable(s)
  be <- matrix(para, ncol = d)   ## be is the matrix of the betas
  be <- cbind(0, be)
  mu1 <- exp(x1 %*% be)
  mu <- mu1 / Rfast::rowsums(mu1) ## fitted values
  ## next we find the fitted values for the compositional vectors with zeros
  ly3 <- ly2
  ind <- which(is.infinite(ly2))
  ly3[ind] <- 0
  mu2 <- exp(x2 %*% be )
  mu2[ind] <- 0
  mu2 <- mu2 / Rfast::rowsums(mu2)
  zeros <-  - sum( lgamma(phi * mu2[mu2>0]), na.rm = TRUE ) + sum( (mu2 * phi - 1) * ly3, na.rm = TRUE )
  ba <- phi * mu
  f <-  - n * lgamma(phi) + sum( lgamma( ba[ba>0] ), na.rm = TRUE ) - sum( (ba - 1) * ly1, na.rm = TRUE ) - zeros
  f
}



