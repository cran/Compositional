compknn.tune <- function(x, ina, nfolds = 10, k = 2:5, type = "S", mesos = TRUE, a = seq(-1, 1, by = 0.1),
                         apostasi = "ESOV", folds = NULL, stratified = TRUE, seed = FALSE, graph = FALSE) {
  n <- dim(x)[1]  ## sample size
  ina <- as.numeric(ina)
  ## number  of nearest neighbours to use
  if ( min(x) == 0 )  a <- a[ a > 0 ]
  if ( is.null(folds) )  folds <- Compositional::makefolds(ina, nfolds = nfolds,
                                                           stratified = stratified, seed = seed)
  nfolds <- length(folds)
  ## The algorithm is repated R times and each time the estimated
  ## percentages are stored in the array per.
  if (apostasi == "ESOV" | apostasi == "taxicab" | apostasi == "CS") {
    a <- a[ a != 0 ]
  }

  if (apostasi == "Ait" | apostasi == "Hellinger" | apostasi == "angular" ) {

    runtime <- proc.time()
    per <- matrix( nrow = nfolds, ncol = length(k) )

    for (vim in 1:nfolds) {
      id <- ina[ folds[[ vim ]] ]  ## groups of test sample
      ina2 <- ina[ -folds[[ vim ]] ]   ## groups of training sample
      g <- Compositional::comp.knn(z[folds[[ vim ]], , drop = FALSE], z[-folds[[ vim ]], ], ina2,
                                     a = NULL, k = k, type = "S", apostasi = apostasi, mesos = mesos)
      be <- g - id
      per[vim, ] <- Rfast::colmeans(be == 0)
    }
    runtime <- proc.time() - runtime

    ela <- Rfast::colmeans(per)
    performance <- max(ela)
    names(performance) <- "rate"
    names(ela) <- paste("k=", k, sep = "")
    best_k <- which.max(ela) + 1

    if (graph)  plot(k, ela, type = "b", xlab = "k nearest neighbours", col = "green", pch = 16,
                     ylab = "Estimated percentage of correct classification", cex.lab = 1.2, cex.axis = 1.2)
    abline(v = k, col = "lightgrey", lty = 2)
    abline(h = seq(min(ela), max(ela), length = 10), col = "lightgrey", lty = 2)

    results <- list(per = ela, performance = performance, best_k = which.max(ela) + 1, runtime = runtime)

  } else {

    per <- array( dim = c(nfolds, length(k), length(a)) )
    runtime <- proc.time()

    for (vim in 1:nfolds) {
      id <- ina[ folds[[ vim ]] ]  ## groups of test sample
      ina2 <- ina[ -folds[[ vim ]] ]   ## groups of training sample
      for ( i in 1:length(a) ) {
        z <- x^a[i]
        z <- x / Rfast::rowsums( z )
        g <- Compositional::comp.knn(z[folds[[ vim ]], , drop = FALSE], z[-folds[[ vim ]], ], ina2,
                                     a = NULL, k = k, type = "S", apostasi = apostasi, mesos = mesos)
        be <- g - id
        per[vim, , i] <- Rfast::colmeans(be == 0)
      }
    }
    runtime <- proc.time() - runtime

    ela <- t( colMeans(per) )
    ## The ela matrix contains the averages of the R
    ## repetitions over alpha and k
    colnames(ela) <- paste("k=", k, sep = "")
    rownames(ela) <- paste("alpha=", a, sep = "")
    ## The code for the heat plot of the estimated percentages
    if (graph)  filled.contour(a, k, ela, ylab = "k nearest-neighbours",
                       xlab = expression(paste(alpha, " values")), cex.lab = 1.2, cex.axis = 1.2 )

    performance <- max(ela)
    names(performance) <- c( "rate")
    confa <- which(ela == performance, arr.ind = TRUE)[1, ]
    results <- list( per = ela, performance = performance, best_a = a[ confa[1] ], best_k = confa[2] + 1,
                     runtime = runtime )
  }
  results
}
