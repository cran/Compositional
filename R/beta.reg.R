beta.reg <- function(y, x, xnew = NULL) {

  regbeta <- function(pa, ly, sly1, x, n) {
    phi <- exp(pa[1])    ;    b <- pa[-1]
    m <- exp( - x %*% b )
    m <- 1 / ( 1 + m )
    a1 <- m * phi   ;   a2 <- phi - a1
    - n * lgamma(phi) + sum( lgamma(a1) ) + sum( lgamma(a2) ) - sum(a1 * ly) - phi * sly1
  }

  x <- model.matrix(y ~ ., data.frame(x) )
  n <- dim(x)[1]
  iniphi <- log( sum( y * (1 - y) ) / Rfast::Var(y) / n )
  ly1 <- log(1 - y)     ;    ly <- log(y) - ly1
  sly1 <- sum(ly1)      ;    sly2 <- sum( log(y) ) + sly1
  suppressWarnings({
    mod1 <- nlm(regbeta, c( iniphi, numeric(dim(x)[2]) ), ly = ly, sly1 = sly1, x = x, n = n, iterlim = 10000 )
    mod2 <- optim(mod1$estimate, regbeta, ly = ly, sly1 = sly1, x = x, n = n, control = list(maxiters = 1000), hessian = TRUE )
  })
  be <- mod2$par[-1]
  se <- sqrt( diag( solve(mod2$hessian) )[-1] )
  stat <- (be/se)^2
  pvalue <- pchisq(stat, 1, lower.tail = FALSE)
  info <- cbind(be, se, stat, pvalue)
  rownames(info) <- colnames(x)

  if ( !is.null(xnew) ) {
    xnew <- model.matrix(~., data.frame(xnew) )
    est <- exp( - as.vector( xnew %*% be[, 1] ) )
    est <- 1 / (1 + est)
  } else  est <- NULL

  list(phi = exp(mod2$par[1]), info = info, loglik = - mod2$value - sly2, est = est)
}
