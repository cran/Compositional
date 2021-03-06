################################
#### BIC for normal mixture models for compositional data
#### Tsagris Michail 5/2015
#### mtsagris@yahoo.gr
#### References: Ryan P. Browne, Aisha ElSherbiny and
#### Paul D. McNicholas (2015)
#### R package mixture: Mixture Models for Clustering and Classification
################################
bic.mixcompnorm <- function(x, G, type = "alr", veo = FALSE, graph = TRUE) {
  ## x is the compositional data
  ## A is the maximum number of components to be considered
  ## type is either 'alr' or 'ilr'
  p <- dim(x)[2]  ## dimensionality of the data
  if (type == "ilr") {
    y0 <- log(x)
    y1 <- y0 - Rfast::rowmeans( y0 )
    y <- tcrossprod( y1, helm(p) )
  } else  y <- log(x[, -1] / x[, 1])

  mod <- mixture::gpcm(y, G = G, mnames = NULL, start = 0, mmax = 100)
  mbic <- mod$BIC[, , 3]  ## BIC for all models
  ## Next, we plot the BIC for all models
  if ( graph ) {
    plot( G, mbic[, 1], type = "b", pch = 9, xlab = "Number of components",
    ylab = "BIC values", ylim = c( min(mbic, na.rm = TRUE), max(mbic, na.rm = TRUE) ),
    cex.lab = 1.3, )
    for ( i in 2:ncol(mbic) )  lines(G, mbic[, i], type = "b", pch = 9, col = i)
  }
  list(mod = mod, BIC = mbic)
}
