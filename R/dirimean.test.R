dirimean.test <- function(x, a) {
  ## x is the compositional data
  ## a is the hypothesized compositional mean vector
  n <- dim(x)[1]  ## sample size
  d <- dim(x)[2] - 1

  if ( min(x) <= 0 || min(a) <= 0 ) {
    res <- paste("There are zeros in the data")

  } else {
    z <- t( log(x) )
    ## loglik is for the 'mle' type
    loglik <- function(phi) {
      phi <- exp(phi)
      ma <- phi * a
      n * lgamma( phi ) - n * sum( lgamma(ma) ) + sum( z * (ma - 1) )
    }

    ## phi under Ho
    phi <- sum(a)
    if ( phi == 1 ) {
      mod0 <- optimize(loglik, c(-20, 20), maximum = TRUE )
      ell0 <- mod0$objective
      phi0 <- exp( mod0$maximum )
      par0 <- phi0 * a

    } else if ( phi > 1 ) {
      ell0 <-  n * lgamma( phi ) - n * sum( lgamma(a) ) + sum( z * ( a - 1 )  )
      par0 <- a
    }
    ## parameters under H1
    mod1 <- diri.nr(x)
    ell1 <- mod1$loglik
    ## test statistic and p-value
    test <- 2 * (ell1 - ell0)
    pvalue <- pchisq(test, d, lower.tail = FALSE)
    param <- rbind(par0, mod1$param)
    rownames(param) <- c("Null", "Alternative")
    if ( is.null( colnames(x) ) ) {
      colnames(param) <- paste("X", 1:c( d + 1 ), sep = " " )
    } else  colnames(param) <- paste("X", 1:c( d + 1 ), sep = " " )

    lik <- c(ell0, ell1)
    names(lik) <- c("Null loglik", "Alternative loglik")
    info <- c(test, pvalue)
    names(info) <- c("Test", "p-value")
    res <- list(param = param, loglik = lik, info = info)

  }

  res

}
