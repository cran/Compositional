\name{rmvt}
\alias{rmvt}
\title{
Multivariate t random values simulation
}

\description{
Multivariate t random values simulation.
}

\usage{
rmvt(n, mu, sigma, v)
}
\arguments{
  \item{n}{
  The sample size, a numerical value.
}
  \item{mu}{
  The mean vector in \eqn{R^d}.
}
  \item{sigma}{
  The covariance matrix in \eqn{R^d}.
}
  \item{v}{
  The degrees of freedom.
}
}

\details{
The algorithm uses univariate normal and chi-square random values and then generates multivariate t values via a spectral decomposition.This function is used by \code{\link{rcompt}}.
}

\value{
A matrix with the simulated data.
}

\references{
Aitchison J. (1986). The statistical analysis of compositional data. Chapman & Hall.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{rdiri}, \link{rmvnorm}, \link{multivt}
}
}

\examples{
x <- as.matrix(iris[, 1:2])
m <- colMeans(x)
s <- var(x)
y <- rmvt(100, m, s, 10)
multivt(y)
}

\keyword{ multivariate t distribution }
\keyword{ random values simulation }
