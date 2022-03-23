\name{Density of compositional data from Gaussian mixture models}
\alias{dmix.compnorm}
\title{
Simulation of compositional data from Gaussian mixture models
}

\description{
Simulation of compositional data from Gaussian mixture models.
}

\usage{
dmix.compnorm(x, mu, sigma, prob, type = "alr", logged = TRUE)
}

\arguments{
\item{x}{
A vector or a matrix with compositional data.
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
The type of trasformation used, either the additive log-ratio ("alr"),
the isometric log-ratio ("ilr") or the pivot coordinate ("pivot") transformation.
}
\item{logged}{
A boolean variable specifying whether the logarithm of the density values to
be returned. It is set to TRUE by default.
}
}

\details{
A sample from a multivariate Gaussian mixture model is generated.
}

\value{
A vector with the density values.

}

\references{
Ryan P. Browne, Aisha ElSherbiny and Paul D. McNicholas (2015). R package mixture: Mixture
Models for Clustering and Classification.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
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
mu <- rowsum(x, ina) / 50
s[, , 1] <- cov(x[ina == 1, ])
s[, , 2] <- cov(x[ina == 2, ])
s[, , 3] <- cov(x[ina == 3, ])
y <- rmixcomp(100, p, mu, s, type = "alr")$x
mod <- dmix.compnorm(y, mu, s, p)
}

\keyword{ random vectors simulation }
\keyword{ Gaussian mixture models }
