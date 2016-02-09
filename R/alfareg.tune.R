################################
#### Tuning the alfa in alfa-regression via K-fold cross validation
#### The bias corrected performance is returned using the
#### Tibshirani and Tibshirani method
#### Tsagris Michail 11/2015
#### mtsagris@yahoo.gr
#### References: Tsagris Michail (2015)
#### Regression analysis with compositional data containing zero values
#### Chilean Journal of Statistics, 6(2): 47-57
#### Tibshirani and Tibshirani (2009),
#### A bias correction for the minimum error rate in cross-validation
#### The Annals of Applied Statistics, 3(1):822-829
################################

alfareg.tune <- function(y, x, a = seq(0.1, 1, by = 0.1), K = 10, nc = 2) {
  ## y is the compositional data (dependent variable)
  ## x is the independent variables
  ## a is a range of values of alpha
  ## K is the number of folds for the K-fold cross validation
  ## nc is how many cores you want to use, default value is 2
  if ( min(y) == 0 )  a <- a[a>0]
  la <- length(a)
  n <- nrow(y)
  nu <- sample(1:n, min( n, round(n / K) * K ) )
  x <- as.matrix(x)
  options(warn = -1)
  mat <- matrix(nu, ncol = K ) # if the length of nu does not fit to a matrix
  ## a warning message should appear
  if (nc == 1) {
    kula <- matrix(nrow = K, ncol = la)
    for (j in 1:la) {
      for (i in 1:K) {
        xu <- x[ mat[, i], ]
        yu <- y[ mat[, i], ]
        xa <- x[ -mat[, i], ]
        ya <- y[ -mat[, i], ]
        mod <- alfa.reg(ya, xa, a[j], xnew = xu)
        yest <- mod$est
        kula[i, j] <- 2 * sum(yu * log(yu / yest), na.rm = T)
      }
    }
    kl <- colMeans(kula)
    opt <- a[ which.min(kl) ]
    val <- which.min(kl)
    per <- min(kl)
    pera <- apply(kula, 1, min)
    bias <- mean( kula[, val] - pera )
  } else {
    options(warn = -1)
    val <- matrix(a, ncol = nc) ## if the length of a is not equal to the
    ## dimensions of the matrix val a warning message should appear
    ## but with options(warn = -1) you will not see it
    cl <- makePSOCKcluster(nc)
    registerDoParallel(cl)
    kula <- foreach(j = 1:nc, .combine = cbind,
    .export = c("alfa.reg", "alfa", "helm", "comp.reg", "multivreg") ) %dopar% {
      ba <- val[, j]
      ww <- matrix(nrow = K, ncol = length(ba) )
      for ( l in 1:length(ba) ) {
        for (i in 1:K) {
          xu <- x[ mat[, i], ]
          yu <- y[ mat[, i], ]
          xa <- x[ -mat[, i], ]
          ya <- y[ -mat[, i], ]
          mod <- alfa.reg(ya, xa, ba[l], xnew = xu)
          yest <- mod$fitted
          ww[i, l] <- 2 * sum(yu * log(yu / yest), na.rm = T)
        }
      }
      return(ww)
    }
    stopCluster(cl)
    kula <- kula[, 1:la]
    kl <- colMeans(kula)
    opt <- a[ which.min(kl) ]
    val <- which.min(kl)
    per <- min(kl)
    pera <- apply(kula, 1, min)
    bias <- mean( kula[, val] - pera )
  }
  plot(a, kula[1, ], type = 'l', ylim = c( min(kula), max(kula) ), xlab = expression(alpha),
  ylab = 'Twice the Kullback Leibler divergence')
  for (i in 2:K)  lines(a, kula[i, ])
  lines(a, kl, col = 2, lty = 2, lwd = 2)
  list(kl = kl, opt = opt, value = per + bias, bias = bias)
}
