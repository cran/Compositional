aitknn.tune <- function (x, ina, nfolds = 10, k = 2:5, mesos = TRUE,
        a = seq(-1, 1, by = 0.1), apostasi = "euclidean", rann = FALSE,
       folds = NULL, stratified = TRUE, seed = NULL, graph = FALSE) {

  if (min(x) == 0) a <- a[a > 0]
  n <- dim(x)[1]
  ina <- as.numeric(ina)
  if ( is.null(folds) )  folds <- Compositional::makefolds(ina, nfolds = nfolds,
                                                           stratified = stratified, seed = seed)
  nfolds <- length(folds)

  runtime <- proc.time()
  if (!rann) {
    ela <- matrix(nrow = length(a), ncol = length(k))
    for (i in 1:length(a)) {
      z <- Compositional::ait(x, a[i], h = FALSE)
      ela[i, ] <- Rfast::knn.cv(folds = folds, nfolds = nfolds, y = ina, x = z,
                  k = k, dist.type = apostasi, type = "C", freq.option = 1)$crit
    }
  } else {
    per <- array(dim = c(nfolds, length(k), length(a)))
    for (i in 1:length(a)) {
      z <- Compositional::ait(x, a[i], h = FALSE)
      for (vim in 1:nfolds) {
        id <- ina[folds[[vim]]]
        ina2 <- ina[-folds[[vim]]]
        aba <- as.vector(folds[[vim]])
        aba <- aba[aba > 0]
        g <- Compositional::ait.knn(z[aba, ], z[-aba, ], ina = ina2, a = NULL, k = k, rann = TRUE)
        be <- g - id
        per[vim, , i] <- Rfast::colmeans(be == 0)
      }
   }
   ela <- t( colMeans(per) )
  }

  colnames(ela) <- paste("k=", k, sep = "")
  rownames(ela) <- paste("alpha=", a, sep = "")
  runtime <- proc.time() - runtime
  if ( graph )  filled.contour(a, k, ela, ylab = "k nearest-neighbours", cex.lab = 1.2, cex.axis = 1.2,
                               xlab = expression(paste(alpha, " values")) )
  opt <- max(ela)
  confa <- as.vector(which(ela == opt, arr.ind = TRUE)[1, ])
  list(ela = ela, performance = max(ela), best_a = a[confa[1]], best_k = confa[2] + 1, runtime = runtime)
}
