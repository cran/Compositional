makefolds <- function(ina, nfolds = 10, stratified = TRUE, seed = FALSE) {
  names <- paste("Fold", 1:nfolds)
  runs <- sapply(names, function(x) NULL)
  if (seed)  set.seed(1234)

  if ( !stratified ) {
    oop <- options(warn = -1) 
	on.exit( options(oop) )
    mat <- matrix( sample(length(ina)), ncol = nfolds )
    for (i in 1:c(nfolds - 1)) runs[[i]] <- mat[, i]
    names <- prod(dim(mat)) - length(ina)
    runs[[nfolds]] <- mat[1:c(nrow(mat) - names), nfolds]
  } else {
    labs <- unique(ina)
    run <- list()
    for (i in 1:length(labs)) {
      names <- which(ina == labs[i])
      run[[i]] <- sample(names)
    }
    run <- unlist(run)
    for (i in 1:length(ina)) {
      k <- i %% nfolds
      if (k == 0)  k <- nfolds
      runs[[k]] <- c(runs[[k]], run[i])
    }
  }
  runs
}
