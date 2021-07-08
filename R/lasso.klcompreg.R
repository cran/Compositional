lasso.klcompreg <- function(y, x, alpha = 1, lambda = NULL, nlambda = 100, type = "grouped", xnew = NULL) {
  mod <- glmnet::glmnet(x, y, alpha = alpha, nlambda = nlambda, lambda = lambda, family = "multinomial",
                        type.multinomial = type)
  est <- NULL
  if ( !is.null(xnew) ) {
    est <- predict(mod, newx = xnew)
  }
  list(mod = mod, est = est)
}
