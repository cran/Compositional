esov <- function(x) {
  Rfast::Dist(x, method = "jensen_shannon")
}
