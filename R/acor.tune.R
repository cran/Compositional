acor.tune <- function(y, x, a = c(-1, 1), type = "dcor") {

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
