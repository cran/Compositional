\name{rmixcomp}
\alias{rmixcomp}
\title{
Simulation of compositional data from Gaussian mixture models
}

\description{
Simulation of compositional data from Gaussian mixture models.
}

\usage{
rmixcomp(n, prob, mu, sigma, type = "alr")
}

\arguments{
  \item{n}{
  The sample size
}
  \item{prob}{
  A vector with mixing  probabilities. Its length is equal to the number of clusters.
}
  \item{mu}{
  A matrix where each row corresponds to the mean vector of each cluster.
}
  \item{sigma}{
  An array consisting of the covariance matrix of each cluster.
}
  \item{type}{
  Should the additive ("type=alr") or the isometric (type="ilr") log-ration be used? The default value is for the additive log-ratio transformation.
}
}

\details{
A sample from a multivariate Gaussian mixture model is generated.
}

\value{
A list including:
  \item{id}{
  A numeric variable indicating the cluster of simulated vector.
}
 \item{x}{
  A matrix containing the simulated compositional data. The number of dimensions will be + 1.
}
}

\references{
Ryan P. Browne, Aisha ElSherbiny and Paul D. McNicholas (2015). R package mixture: Mixture Models for Clustering and Classification.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{mix.compnorm}, \link{bic.mixcompnorm}
}
}

\examples{
p <- c(1/3, 1/3, 1/3)
mu <- matrix(nrow = 3, ncol = 4)
s <- array( dim = c(4, 4, 3) )
x <- as.matrix(iris[, 1:4])
ina <- as.numeric(iris[, 5])
for (i in 1:3) {
  mu[i, ] <- colMeans(x[ina == i, ])
  s[, , i] <- cov(x[ina == i, ])
}
y <- rmixcomp(100, p, mu, s, type = "alr")
}

\keyword{ random vectors simulation }
\keyword{ Gaussian mixture models }
