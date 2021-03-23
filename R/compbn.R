compbn <- function(x, type = "fedhc", max_k = 3, alpha = 0.05,
                   robust = FALSE, ini.stat = NULL, R = NULL,
                   restart = 10, tabu = 10, score = "bic-g", blacklist = NULL, whitelist = NULL) {

  if ( type == "fedhc" ) {
    mod <- pchc::fedhc(x, alpha = alpha, robust = robust, ini.stat = ini.stat, R = R,
                       restart = restart, score = score, blacklist = blacklist, whitelist = whitelist)

  } else if ( type == "pchc" ) {
    mod <- pchc::pchc(x, alpha = alpha, robust = robust, ini.stat = ini.stat, R = R,
                       restart = restart, score = score, blacklist = blacklist, whitelist = whitelist)

  } else if ( type == "mmhc" ) {
    mod <- pchc::mmhc(x, max_k = max_k, alpha = alpha, robust = robust, ini.stat = ini.stat, R = R,
                       restart = restart, score = score, blacklist = blacklist, whitelist = whitelist)

  } else if ( type == "fedtabu" ) {
    mod <- pchc::fedtabu(x, alpha = alpha, robust = robust, ini.stat = ini.stat,
                       R = R, tabu = tabu, score = score, blacklist = blacklist, whitelist = whitelist)

  } else if ( type == "pctabu" ) {
    mod <- pchc::pctabu(x, alpha = alpha, robust = robust, ini.stat = ini.stat,
           R = R, tabu = tabu, score = score, blacklist = blacklist, whitelist = whitelist)

  } else if ( type == "mmtabu" ) {
    mod<- pchc::mmtabu(x, max_k = 3, alpha = alpha, robust = robust, ini.stat = ini.stat,
                       R = R, tabu = tabu, score = score, blacklist = blacklist, whitelist = whitelist)
  }

  mod
}
