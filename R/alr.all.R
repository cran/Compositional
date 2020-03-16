alr.all <- function(x) {
  y1 <- Rfast::Log(x)
  d <- dim(x)[2]
  y <- y1[, -1] - y1[, 1]
  for (i in 2:d)  y <- cbind(y, y1[, -c(1:i)] - y1[, i])
  y
}
