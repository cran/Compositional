\name{Multivariate t random values simulation on the simplex}
\alias{rcompt}
\title{
Multivariate t random values simulation on the simplex
}

\description{
Multivariate t random values simulation on the simplex.
}

\usage{
rcompt(n, m, s, dof, type = "alr")
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
\item{dof}{
The degrees of freedom.
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
Aitchison J. (1986). The statistical analysis of compositional data. Chapman & Hall.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} and Giorgos Athineou <gioathineou@gmail.com>
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
m <- Rfast::colmeans(x)
s <- var(x)
y <- rcompt(100, m, s, 10)
comp.den(y, dist = "t")
ternary(y)
}

\keyword{ multivariate t distribution }
\keyword{ random values simulation }
