alfalasso.tune <- function(y, x, a = seq(-1, 1, by = 0.1), model = "gaussian", lambda = NULL,
                           type.measure = "mse", nfolds = 10, folds = NULL, stratified = FALSE) {

  n <- length(y)  ;   f <- numeric(n)
  if ( is.null(folds) )  folds <- Compositional::makefolds(y, nfolds = nfolds, stratified = stratified)
  for ( i in 1:nfolds )   f[ folds[[ i ]] ] <- i
  nfolds <- length(folds)

  res <- matrix(nrow = length(a), ncol = 2)
  rownames(res) <- paste("alpha=", a, sep = "")
  colnames(res) <- c("lambda", "performance")

  for ( i in 1:length(a) ) {
    xa <- Compositional::alfa(x, a[i], h = TRUE)$aff ## apply the alpha-transformation
    mod <- glmnet::cv.glmnet(xa, y, family = model, lambda = lambda, type.measure = type.measure,
                             nfolds = nfolds, foldid = f)
    res[i, ] <- c( mod$lambda[ mod$index[1] ], mod$cvm[ mod$index[1] ])
  }
  res
}
