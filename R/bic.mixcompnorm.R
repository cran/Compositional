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

  if ( type == "alr" ) {
    y <- log(x[, -1] / x[, 1])
  } else if (type == "ilr") {
    y0 <- log(x)
    y1 <- y0 - Rfast::rowmeans( y0 )
    y <- tcrossprod( y1, Compositional::helm(p) )
  } else if ( type == "pivot" ) {
    y <- Compositional::pivot(x)
  }

  mod <- mixture::gpcm(y, G = G, mnames = NULL, start = 0, mmax = 100, veo = veo)
  mbic <- mod$BIC[, , 3, drop = FALSE]  ## BIC for all models
  ## Next, we plot the BIC for all models
  if ( graph ) {
    plot( G, mbic[, 1, drop = FALSE], type = "b", pch = 9, xlab = "Number of components",
    ylab = "BIC values", ylim = c( min(mbic, na.rm = TRUE), max(mbic, na.rm = TRUE) ),
    cex.lab = 1.2, cex.axis = 1.2, xaxt = "n")
    axis(1, at = G, labels = G)
    abline(v = G, col = "lightgrey", lty = 2)
  	abline(h = seq( min(mbic), max(mbic), length = 10 ), col = "lightgrey", lty = 2)
    for ( i in 2:ncol(mbic) )  lines(G, mbic[, i], type = "b", pch = 9, col = i)
  }

  pou <- which( mbic == max(mbic, na.rm = TRUE),arr.ind = TRUE )
  if ( length(G) == 1 ) {
    optG <- G
  } else  optG <- G[ pou[1, 1] ]
  optmodel <- colnames(mbic)[ pou[1, 2] ]

  list(mod = mod, BIC = mbic, optG = optG, optmodel = optmodel)
}
