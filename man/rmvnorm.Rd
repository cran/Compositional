\name{rmvnorm}
\alias{rmvnorm}
\title{
Multivariate normal random values simulation
}

\description{
Multivariate normal random values simulation.
}

\usage{
rmvnorm(n, mu, sigma)
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
}

\details{
The algorithm uses univariate normal random values and transforms them to multivariate via a spectral decomposition.
It is slower than the command \code{\link{mvrnorm}} available from MASS, but it allows for singular covariance matrices.
This function is used by \code{\link{rcompnorm}}.
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
\code{\link{comp.den}, \link{rdiri}, \link{rmvt}
}
}

\examples{
x <- as.matrix(iris[, 1:2])
m <- colMeans(x)
s <- var(x)
y <- rmvnorm(100, m, s)
colMeans(y)
var(y)
}

\keyword{ multivariate normal distribution }
\keyword{ random values simulation }
