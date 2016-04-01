\name{rcompnorm}
\alias{rcompnorm}
\title{
Multivariate normal random values simulation on the simplex
}

\description{
Multivariate normal random values simulation on the simplex.
}

\usage{
rcompnorm(n, m, s, type = "alr")
}
\arguments{
  \item{n}{
  The sample size, a numerical value.
}
  \item{m}{
  The mean vector in \eqn{R^d}.
}
  \item{s}{
  The covariance matrix in \eqn{R^d}.
}
  \item{type}{
  The alr (type = "alr") or the ilr (type = "ilr") is to be used for closing the Euclidean data onto the simplex.
}
}

\details{
The algorithm is straightforward, generate random values from a multivariate normal distribution in \eqn{R^d} and brings the
values to the simplex \eqn{S^d} using the inverse of a log-ratio transformation.
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
\code{\link{comp.den}, \link{rdiri}, \link{rcompt}, \link{rcompsn}
}
}

\examples{
x <- as.matrix(iris[, 1:2])
m <- colMeans(x)
s <- var(x)
y <- rcompnorm(100, m, s)
comp.den(y)
ternary(y)
}

\keyword{ multivariate normal distribution }
\keyword{ random values simulation }
