lasso.compreg <- function(y, x, alpha = 1, lambda = NULL, nlambda = 100, xnew = NULL) {
  y <- Compositional::alr(y)
  mod <- glmnet::glmnet(x, y, alpha = alpha, nlambda = nlambda, lambda = lambda, family = "mgaussian")
  est <- NULL
  if ( !is.null(xnew) ) {
    est <- predict(mod, newx = xnew)
  }
  list(mod = mod, est = est)
}
