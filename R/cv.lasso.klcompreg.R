cv.lasso.klcompreg <- function(y, x, alpha = 1, type = "grouped", nfolds = 10,
                               folds = NULL, seed = NULL, graph = FALSE) {
  n <- dim(y)[1]  ## sample size
  ina <- 1:n
  if ( is.null(folds) )  folds <- Compositional::makefolds(ina, nfolds = nfolds,
                                                           stratified = FALSE, seed = seed)
  nfolds <- length(folds)
  foldid <- numeric(n)
  for ( i in 1:nfolds )  foldid[ folds[[ i ]] ] <- i

  mod <- glmnet::cv.glmnet(x, y, alpha = alpha, family = "multinomial", type.multinomial = type,
                           foldid = foldid, type.measure = "deviance")
  if ( graph )  plot(mod, cex.lab = 1.2, cex.axis = 1.2)
  mod
}
