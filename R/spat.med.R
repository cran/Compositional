################################
#### Spatial median
#### Tsagris Michail 10/2014
#### References: Jyrki Mottonen, Klaus Nordhausen and Hannu Oja (2010)
#### Asymptotic theory of the spatial median
#### In Nonparametrics and Robustness in Modern Statistical Inference and Time Series
#### Analysis: A Festschrift in honor of Professor Jana Jureckova
#### On computation of spatial median for robust data mining (2005)
#### T. Karkkaminen and S. Ayramo
#### Evolutionary and Deterministic Methods for Design, Optimization
#### and Control with Applications to Industrial and Societal Problems
#### EUROGEN 2005
#### R. Schilling, W.Haase, J. Periaux, H. Baier, G. Bugeda (Eds)
#### FLM, Munich, 2005
#### http://users.jyu.fi/~samiayr/pdf/ayramo_eurogen05.pdf
#### mtsagris@yahoo.gr
################################

spat.med <- function(x) {
  ## contains the data

  x <- as.matrix(x)
  y <- t(x)
  u1 <- as.vector( Rfast::colMedians(x) )

  z <- y - u1
  ww <- 1 / sqrt( colSums(z^2) )

  wei <- ww / sum(ww)
  u2 <-  as.vector( y %*% wei )

  while ( sum( abs(u2 - u1) ) > 1e-9 ) {
    z <- y - u2
    u1 <- u2

    ww <- 1 / sqrt( colSums(z^2) )
    wei <- ww / sum(ww)
    u2 <- as.vector( y %*% wei )

    if ( any( is.na( u2 ) ) ) {
      u2 <- u1
    }

  }

  u2

}
