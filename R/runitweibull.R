runitweibull <- function(n, a, b) {
  exp( -rweibull(n, a, b) )
}
