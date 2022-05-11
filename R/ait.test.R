ait.test <- function(x1, x2, type = 1, alpha = 0.05) {

  y1 <- Compositional::alr(x1)
  y2 <- Compositional::alr(x2)
  m1 <- Rfast::colmeans(y1)
  m2 <- Rfast::colmeans(y2)
  d <- dim(y1)[2]
  n1 <- dim(y1)[1]  ;  n2 <- dim(y2)[1]
  s1 <- ( crossprod(y1) - n1 * tcrossprod(m1) ) / n1
  s2 <- ( crossprod(y2) - n2 * tcrossprod(m2) ) / n2

  if ( type == 1 ) {
    sp <- (n1 * s1 + n2 * s2) / (n1 + n2)
    sc <- sp + n1 * n2 / (n1 + n2) * sum( (m1 - m2)^2 )
    dof <- 0.5 * d * (d + 3)
    stat <- n1 * log( det(sc) / det(s1) ) + n2 * log( det(sc) / det(s2) )
    pvalue <- pchisq(stat, dof, lower.tail = FALSE)

  } else if ( type == 2 ) {
    sp <- (n1 * s1 + n2 * s2) / (n1 + n2)
    dof <- 0.5 * d * (d + 1)
    stat <- n1 * log( det(sp) / det(s1) ) + n2 * log( det(sp) / det(s2) )
    pvalue <- pchisq(stat, dof, lower.tail = FALSE)

  } else {
    i <- 1
    s1h <- s1  ;  s2h <- s2
    s1inv <- solve(s1h)  ;  s2inv <- solve(s2h)
    mha <- solve( n1 * s1inv + n2 * s2inv, n1 * s1inv %*% m1 + n2 * s2inv %*% m2)
    s1h <- s1h + tcrossprod(m1 - mha)
    s2h <- s2h + tcrossprod(m2 - mha)
    s1inv <- solve(s1h)  ;  s2inv <- solve(s2h)
    mhb <- solve( n1 * s1inv + n2 * s2inv, n1 * s1inv %*% m1 + n2 * s2inv %*% m2 )
    while ( sum( abs(mha - mhb) ) > 1e-6 ) {
      i <- i + 1
      mha <- mhb
      mhb <- solve( n1 * s1inv + n2 * s2inv, n1 * s1inv %*% m1 + n2 * s2inv %*% m2 )
      s2h <- s1h + tcrossprod(m1 - mhb)
      s2h <- s1h + tcrossprod(m2 - mhb)
    }
    dof <- d
    stat <- n1 * log( det(s1h) / det(s1) ) + n2 * log( det(s2h) / det(s2) )
    pvalue <- pchisq(stat, dof, lower.tail = FALSE)
  }

  crit <- qchisq(1 - alpha, dof)
  info <- c(stat, pvalue, crit, dof)
  names(info) <- c("test", "p-value", "critical", "degrees of freedom")
  info

}
