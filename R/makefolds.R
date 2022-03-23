makefolds <- function(ina, nfolds = 10, stratified = TRUE, seed = NULL) {
  names <- paste("Fold", 1:nfolds)
  runs <- sapply(names, function(x) NULL)
  if ( !is.null(seed) )  set.seed(seed)

  if ( !stratified ) {
    oop <- options(warn = -1)
    rat <- length(ina) %% nfolds
    mat <- matrix( sample( length(ina) ), ncol = nfolds )
    mat[-c( 1:length(ina) )] <- NA
    on.exit( options(oop) )
    for ( i in 1:c(nfolds - 1) )  runs[[ i ]] <- mat[, i]
    a <- prod(dim(mat)) - length(ina)
    runs[[ nfolds ]] <- mat[1:c(nrow(mat) - a), nfolds]
  } else {
    labs <- unique(ina)
    run <- list()
    for (i in 1:length(labs)) {
      names <- which( ina == labs[i] )
      run[[ i ]] <- sample(names)
    }
    run <- unlist(run)
    for ( i in 1:length(ina) ) {
      k <- i %% nfolds
      if ( k == 0 )  k <- nfolds
      runs[[ k ]] <- c( runs[[ k ]], run[i] )
    }
  }
  
  for (i in 1:nfolds)  {
    if ( any( is.na(runs[[ i ]]) ) )  runs[[ i ]] <- runs[[ i ]][ !is.na(runs[[ i ]]) ]
  }
  
  if ( length(runs[[ nfolds ]]) == 0 ) runs[[ nfolds ]] <- NULL
  runs
}


