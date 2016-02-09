\name{rcompsn}
\alias{rcompsn}
\title{
Multivariate skew normal random values simulation on the simplex
}
\description{
Multivariate skew normal random values simulation on the simplex.
}

\usage{
rcompsn(n, xi, Omega, alpha, dp = NULL, type = "alr")
}

\arguments{
  \item{n}{
  The sample size, a numerical value.
}
  \item{xi}{
  A numeric vector of length \eqn{d} representing the location parameter of the distribution.
}
  \item{Omega}{
  A \eqn{d \times d} symmetric positive-definite matrix of dimension.
}
  \item{alpha}{
  A numeric vector which regulates the slant of the density.
}
  \item{dp}{
  A list with three elements, corresponding to xi, Omega and alpha described above. The default value is FALSE.
  If dp is assigned, individual parameters must not be specified.
}
  \item{type}{
  The alr (type = "alr") or the ilr (type = "ilr") is to be used for closing the Euclidean data onto the simplex.
}
}

\details{
The algorithm is straightforward, generate random values from a multivariate t distribution in \eqn{R^d} and brings the
values to the simplex \eqn{S^d} using the inverse of a log-ratio transformation.
}

\value{
A matrix with the simulated data.
}
\references{
Azzalini, A. and Dalla Valle, A. (1996). The multivariate skew-normal distribution. Biometrika, 83(4): 715-726.

Azzalini, A. and Capitanio, A. (1999). Statistical applications of the multivariate skew normal distribution. J.Roy.Statist.Soc. B, 61(3):579-602. Full-length version available at http://arXiv.org/abs/0911.2093

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
\code{\link{comp.den}, \link{rdiri}, \link{rcompnorm}, \link{rmvt}
}
}

\examples{
x <- as.matrix(iris[, 1:2])
par <- sn::msn.mle(y = x)$dp
y <- rcompsn(100, dp = par)
comp.den(y, dist = "skewnorm")
ternary(y)
}

\keyword{ multivariate skew normal distribution }
\keyword{ random values simulation }
