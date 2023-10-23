lcglm.aov <- function(mod0, mod1) {
  stat <- mod0$devi - mod1$devi
  pvalue <- pchisq(stat, 1, lower.tail = FALSE)
  res <- c(stat, pvalue)
  names(res) <- c( "Chi-square statistic", "p-value" )
  res
}
