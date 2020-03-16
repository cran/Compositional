comp.test <- function(x, ina, test = "james", R = 0, ncores = 1, graph = FALSE) {
  ## x contains all the groups together
  ## ina is the group indicator variable
  ## test is the type of test to be used
  ## R takes values 0, 1, 2, or much higher, like 999.
  ## If test is "hotel", "maov" or "maovjames" the value of
  ## R is not taken into account
  ## If test is "james", R can be either 1 (James test)
  ## or 2 (MNV modification of the James test).
  ## If R is 0, it becomes 1 by default.
  ## If test is "el" or "eel", R can be either 0, 1 or 2. The value of 0 means
  ## that the asymptotic chi-squre distribution is used. The value of 1
  ## means that the James corrected chi-square distribution is used.
  ## The value of 2 means that the F distribution used in the MNV test
  ## is used.
  ## if R>2 bootstrap calculation of the p-value is performed
  ## 999 bootstrap resamples are set by default
  ## bootstrap is used for the p-value
  ## ncores is the number of cores you want to use
  ## requires(doParallel)
  ## if graph is TRUE, the bootstrap statics are plotted
  y <- alfa(x, 1)$aff  ## the alpha-transformation with alpha = 1
  ina <- as.numeric(ina)  ## the group indicator variable
  k <- max(ina)  ## the number of groups
  ## default value in the case of MANOVA and bad specification of test
  ## is the James MANOVA
  if ( k > 2 & ( test != "maovjames" || test != "maov" ) )   result <- Compositional::maovjames(x, ina)
  ## multi-sample case
  if (k > 2) {

    if (test == "maov") {
      result <- Compositional::maov(x, ina)
    } else if ( test == "maovjames" )  result <- Compositional::maovjames(x, ina)
    ## two sample case
  } else if ( k == 2 ) {

    if ( test == "hotel" ) {
      result <- Compositional::hotel2T2(y[ina == 1, ], y[ina == 2, ], R = R, graph = graph)
    } else if ( test == "james" ) {
      result <- Compositional::james(y[ina == 1, ], y[ina == 2, ], R = R, graph = graph)
    } else if ( test == "el" ) {
      result <- Compositional::el.test2(y[ina == 1, ], y[ina == 2, ], R = R, ncores = ncores, graph = graph)
    } else if ( test == "eel" ) {
      result <- Compositional::eel.test2(y[ina == 1, ], y[ina == 2, ], R = R, graph = graph)
    }
  }
  result
}
