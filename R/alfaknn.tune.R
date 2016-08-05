###############################
#### Classification for compositional data using the alpha-transformation
#### Tuning the k-NN algorithm
#### Tsagris Michail 8/2015
#### References: Tsagris, M., Preston S. and Wood A.T.A. (2016).
#### Improved classication for compositional data using the alpha-transformation
#### Journal of Classification (To appear)
#### http://arxiv.org/pdf/1506.04976v2.pdf
#### mtsagris@yahoo.gr
################################

alfaknn.tune <- function(x, ina, M = 10, A = 5, type = "S", mesos = TRUE,
  a = seq(-1, 1, by = 0.1), mat = NULL, graph = FALSE) {

  ## x is the matrix containing the data
  ## A is the maximum number of neighbours to use
  ## ina indicates the groups, numerical variable
  ## a is a vector containing the values of the power parameter
  ## type is either 'S' or 'NS'. Should the standard k-NN be use or not
  ## if mesos is TRUE, then the arithmetic mean distange of the k nearest
  ## points will be used.
  ## If not, then the harmonic mean will be used. Both of these apply for
  ## the non-standard algorithm, that is when type='NS'

  x <- as.matrix(x)  ## makes sure the x is a matrix
  x <- x / as.vector( Rfast::rowsums(x) )  ## makes sure the the data sum to 1
  if ( min(x) == 0 )  a <- a[a>0]  ## checks for any zeros in the data
  n <- nrow(x)  ## sample size
  if ( A >= min(table(ina)) )  A <- min(table(ina)) - 3  ## The maximum
  ## number of nearest neighbours to use
  ina <- as.numeric(ina) ## makes sure ina is numeric
  ng <- max(ina)  ## The number of groups

  ## as the one specified by the user
  ## The next two functions split the sample into R different test
  ## and training datasets
  ## The test dataset is chosen via stratified or simple random sampling
  ## will be stored in the array called per
  ## if seed==TRUE then the results will always be the same

  dis <- matrix(0, n, n)
  ## The next two functions split the sample into R different test
  ## and training datasets
  ## The test dataset is chosen via stratified or simple random sampling
  ## will be stored in the array called per

  if ( is.null(mat) ) {
    nu <- sample(1:n, min( n, round(n / M) * M ) )
    ## It may be the case this new nu is not exactly the same
    ## as the one specified by the user
    ## to a matrix a warning message should appear
    options(warn = -1)
    mat <- matrix( nu, ncol = M ) # if the length of nu does not fit
  } else  mat <- mat

  M <- ncol(mat)
  rmat <- nrow(mat)

  ## The algorithm is repeated R times and each time the estimated
  ## percentages are stored in the array per.

    runtime <- proc.time()
    per <- array( dim = c( M, A - 1, length(a) ) )  ## The estimated percentages

    for ( i in 1:length(a) ) {
     dis <- alfadist(x, a[i]) ## euclidean distance matrix to the
      ## alpha-transformed data
      ## The k-NN algorith is calculated R times. For every repetition a
      ## test sample is chosen and its observations are classified

      for (vim in 1:M) {

        id <- as.vector( ina[ mat[, vim] ] )  ## groups of test sample
        ina2 <- as.vector( ina[ -mat[, vim] ] )   ## groups of training sample
        aba <- as.vector( mat[, vim] )
        aba <- aba[aba > 0]
        apo <- dis[aba, -aba]
        ta <- matrix(nrow = rmat, ncol = ng)

        if (type == "NS") {
          ## Non Standard algorithm
          for ( j in 1:c(A - 1) ) {
            knn <- j + 1
            for (l in 1:ng) {
              dista <- apo[, ina2 == l]
              dista <- t( apply(dista, 1, sort) )
              if (mesos == TRUE) {
                ta[, l] <- as.vector( Rfast::rowmeans( dista[, 1:knn] ) )
              } else {
                ta[, l] <- knn / as.vector( Rfast::rowsums( 1 / dista[, 1:knn] ) )
              }
            }
            g <- apply(ta, 1, which.min)
            per[vim, j, i] <- mean(g == id)
          }

        } else if (type == "S") {
          ## Standard algorithm
          for (j in 1:c(A - 1) ) {
            g <- numeric(rmat)
            knn <- j + 1
            for (k in 1:rmat) {
              xa <- cbind(ina2, apo[k, ])
              qan <- xa[order(xa[, 2]), ]
              sa <- qan[1:knn, 1]
              tab <- table(sa)
              g[k] <- as.integer(names(tab)[which.max(tab)])
            }
            per[vim, j, i] <- mean(g == id)
          }
        }
      }
    }

    ela <- matrix(nrow = length(a), ncol = A - 1)
    for ( i in 1:length(a) )  ela[i, ] <- colMeans(per[, , i])
    ## The ela matrix contains the averages of the R
    ## repetitions over alpha and k
    colnames(ela) <- paste("k=", 2:A, sep = "")
    rownames(ela) <- paste("alpha=", a, sep = "")

    ## The code for the heat plot of the estimated percentages
    if (graph == TRUE) {
      fields::image.plot(a, 2:A, ela, col = grey(1:11/11),
                         ylab = "k nearest-neighbours",
                         xlab = expression(paste(alpha, " values")) )
    }

    opt <- max(ela)
    confa <- as.vector( which(ela == opt, arr.ind = TRUE)[1, ] )
    bias <- numeric(M)
    for (i in 1:M) {
      bias[i] <- opt - per[ i, confa[2], confa[1] ]
    }
    bias <- mean(bias)
    performance <- c(opt - bias, bias)
    names(performance) <- c( "rate", "bias" )
    runtime <- proc.time() - runtime

    list( ela = ela, performance = performance, best_a = a[ confa[1] ],
    best_k = confa[2] + 1, runtime = runtime )

}
