compknn.tune <- function(x, ina, M = 10, A = 5, type = "S", mesos = TRUE, a = seq(-1, 1, by = 0.1),
                         apostasi = "ESOV", mat = NULL, graph = FALSE) {
  n <- dim(x)[1]  ## sample size
  ina <- as.numeric(ina)
  if ( A >= min(table(ina)) )  A <- min( table(ina) ) - 3  ## The maximum
  ## number  of nearest neighbours to use
  if ( min(x) == 0 )  a <- a[ a > 0 ]
  if ( is.null(mat) ) {
    nu <- sample(1:n, min( n, round(n / M) * M ) )
    ## It may be the case this new nu is not exactly the same
    ## as the one specified by the user
    ## to a matrix a warning message should appear
    options(warn = -1)
    mat <- matrix( nu, ncol = M ) # if the length of nu does not fit
  } else  mat <- mat

  M <- dim(mat)[2]
  ## The algorithm is repated R times and each time the estimated
  ## percentages are stored in the array per.
  if (apostasi == "ESOV" | apostasi == "taxicab" | apostasi == "CS") {
    a <- a[ a != 0 ]
  }
  per <- array( dim = c(M, A - 1, length(a)) )
  runtime <- proc.time()

  for (vim in 1:M) {
    id <- ina[ mat[, vim] ]  ## groups of test sample
    ina2 <- ina[ -mat[, vim] ]   ## groups of training sample
    aba <- as.vector( mat[, vim] )
    aba <- aba[aba > 0]
    for ( i in 1:length(a) ) {
	  z <- x^a[i]
      z <- x / Rfast::rowsums( z )
      g <- comp.knn(z[aba, ], z[-aba, ], ina2, a = NULL, k = 2:A, type = "S", apostasi = apostasi, mesos = mesos)
      be <- g - id
      per[vim, , i] <- Rfast::colmeans(be == 0)
    }
  }
  runtime <- proc.time() - runtime

  if (apostasi == "Ait" | apostasi == "Hellinger" | apostasi == "angular" ) {

    ela <- Rfast::colmeans(per)
    performance <- max(ela)
    names(performance) <- "rate"
    names(ela) <- paste("k=", 2:A, sep = "")
    best_k <- which.max(ela) + 1

    if (graph)  plot(2:A, ela, type = "b", xlab = "k nearest neighbours", pch = 9, col = 2,
                     ylab = "Estimated percentage of correct classification", cex.lab = 1.3)
    results <- list(ela = ela, performance = performance, best_k = which.max(ela) + 1, runtime = runtime)

  } else {
    ela <- t( colMeans(per) )
    ## The ela matrix contains the averages of the R
    ## repetitions over alpha and k
    colnames(ela) <- paste("k=", 2:A, sep = "")
    rownames(ela) <- paste("alpha=", a, sep = "")
    ## The code for the heat plot of the estimated percentages
    if (graph)  fields::image.plot(a, 2:A, ela, col = grey(1:11/11), ylab = "k nearest-neighbours",
                                   xlab = expression(paste(alpha, " values")), cex.lab = 1.3 )

    performance <- max(ela)
    names(performance) <- c( "rate")
    confa <- which(ela == performance, arr.ind = TRUE)[1, ]
    results <- list( ela = ela, performance = performance, best_a = a[ confa[1] ], best_k = confa[2] + 1, runtime = runtime )
  }
  results
}
