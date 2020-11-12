acor.tune <- function(y, x, a = seq(-1, 1, by = 0.1), type = "dcor") {

  if ( min(x) == 0  |  min(y) == 0 )  a = a[ a > 0 ]

  runtime <- proc.time()

  if ( type == "cancor1" ) {
    pa <- function(a, y, x)  sum( log( 1 - acor( y, x, a, type = "cancor" )^2 ) )
    alfa <- optimize( pa, a, y = y, x = x, maximum = TRUE )
  } else if ( type == "cancor2" ) {
    pa <- function(a, y, x)  acor( y, x, a, type = "cancor" )[1]
    alfa <- optimize( pa, a, y = y, x = x, maximum = TRUE )
  } else if ( type == "dcor" ) {
    pa <- function(a, y, x)  acor( y, x, a, type = "dcor" )
    alfa <- optimize( pa, a, y = y, x = x, maximum = TRUE )
  }

  runtime <- proc.time() - runtime

  list( alfa = alfa$maximum, acor = alfa$objective, runtime = runtime )
}
