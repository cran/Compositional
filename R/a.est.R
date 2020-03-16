a.est <- function(x) {
  ## x contains compositional data
   runtime <- proc.time()
   opt <- optimize( a.mle, c(-1, 1), x = x, maximum = TRUE )
   best <- opt$maximum
   mod <- Compositional::alpha.mle(x, best)
   runtime <- proc.time() - runtime
   list(runtime = runtime, best = best, loglik = mod$loglik,
        p = mod$p, mu = mod$mu, su = mod$su)
}
