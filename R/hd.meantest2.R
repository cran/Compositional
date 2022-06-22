hd.meantest2 <- function(y1, y2, R = 1) {

  z1 <- Rfast::Log(y1)
  z1 <- z1 - Rfast::rowmeans( z1 )
  z2 <- Rfast::Log(y2)
  z2 <- z2 - Rfast::rowmeans( z2 )
  n1 <- dim(y1)[1]   ;   n2 <- dim(y2)[1]   ;    p <- dim(y1)[2]
  m1 <- Rfast::colmeans(z1)
  m2 <- Rfast::colmeans(z2)
  s1 <- (n1 - 1) * Rfast::colVars(z1)
  s2 <- (n2 - 1) * Rfast::colVars(z2)
  gii <- (s1 + s2) / (n1 + n2)
  Mn <- (m1 - m2)^2 / gii
  Mn <- (n1 * n2)/(n1 + n2) * max(Mn)
  pvalue <- 1 - exp( -1/sqrt(pi) * exp( ( - Mn + 2 * log(p) - log( log(p) ) )/2 ) )   
  
  if ( R > 1 ) {
    mc <- (n1 * m1 / s1 + n2 * m2 / s2) / (n1/s1 + n2/s2)
    mc1 <-  - m1 + mc	
    mc2 <-  - m2 + mc  
    x1 <- Rfast::eachrow(z1, mc1, oper = "+")
    x2 <- Rfast::eachrow(z2, mc2, oper = "+")
    Mnb <- numeric(R)  
    for (i in 1:R) {
      xb1 <- x1[Rfast2::Sample.int(n1, n1, replace = TRUE), ]   	
      xb2 <- x2[Rfast2::Sample.int(n2, n2, replace = TRUE), ]
      m1 <- Rfast::colmeans(xb1)
      m2 <- Rfast::colmeans(xb2)
      s1 <- (n1 - 1) * Rfast::colVars(xb1)
      s2 <- (n2 - 1) * Rfast::colVars(xb2)
      gii <- (s1 + s2) / (n1 + n2)
      stat <- (m1 - m2)^2 / gii
      Mnb[i] <- max(stat)
    }
    pvalue <- ( sum( (n1 * n2)/(n1 + n2) * Mnb >= Mn ) + 1 ) / (R + 1)
  }

  res <- c(Mn, pvalue)
  names(res) <- c("Mn", "p-value")
  res 
  
} 