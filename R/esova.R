## computes ESOV distance between xnew and x
esova <- function(xnew, x) {
  Rfast::dista(xnew, x, type="jensen_shannon")
}


