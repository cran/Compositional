################################
#### BIC for normal mixture models for compositional data
#### Tsagris Michail 5/2015
#### mtsagris@yahoo.gr
#### References: Ryan P. Browne, Aisha ElSherbiny and
#### Paul D. McNicholas (2015)
#### R package mixture: Mixture Models for Clustering and Classification
################################
bic.alfamixnorm <- function(x, G, a = seq(-1, 1, by = 0.1), veo = FALSE, graph = TRUE) {
  ## x is the compositional data
  ## A is the maximum number of components to be considered
  ## type is either 'alr' or 'ilr'
  n <- dim(x)[1]   ;   p <- dim(x)[2]  ## dimensionality of the data
  if ( min(x) == 0 )  a <- a[a > 0]
  names <- paste("Alfa=", a)
  abic <- sapply(names, function(x) NULL)

  for ( i in 1:length(a) ) {
    z <- alfa(x, a[i])
    y <- z$aff
    sa <- z$sa
    mod <- mixture::gpcm(y, G = G, mnames = NULL, start = 0, mmax = 100, veo = veo)
    abic[[ i ]] <- mod$BIC[, , 3] + 2 * sa - log(n) ## BIC for all models
    pou <- which( is.na(abic[[ i ]] ) )
    if ( length(pou) > 0 )  abic[[ i ]][pou] <-  -Inf
  }
  ## Next, we plot the BIC for all models
  if ( graph ) {
    bica <- matrix(nrow = length(G), ncol = length(a) )
    for ( i in 1:length(a) )  bica[, i] <- Rfast::rowMaxs( abic[[ i ]], value = TRUE )
    plot( a, bica[1, ], type = "b", pch = 9, xlab = expression( paste(alpha, " values") ),
          ylab = "BIC values", ylim = c( min(bica, na.rm = TRUE), max(bica, na.rm = TRUE) ),
          cex.lab = 1.3, cex.axis = 1.2)
    for ( i in 2:max(G) )  lines(a, bica[i, ], type = "b", pch = 9, col = i)
  }
  abic
}
