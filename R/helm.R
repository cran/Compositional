################################
#### The Helmert sub-matrix
#### Tsagris Michail 5/2011  
#### References: John Aitchison (2003) 
#### The Statistical Analysis of Compositional Data p. 99 Blackburn Press 
#### Lancaster H. O. (1965). The Helmert matrices. 
#### The American Mathematical Monthly 72(1): 4-12.
################################

helm <- function(n) {
  h <- matrix(numeric(n^2), nrow = n)
  h[1, ] <- 1/sqrt(n)
  for (i in 2:n) {
    for (j in 1:i - 1) h[i, j] <- 1/sqrt(i * (i - 1))
    h[i, j + 1] <- -sum(h[i, ])
  }
  h[c(2:n), ]
}