\name{Estimation of the value of alpha via the profile log-likelihood}
\alias{alfa.profile}
\title{
Estimation of the value of \eqn{\alpha} via the alfa profile log-likelihood
}

\description{
Estimation of the value of \eqn{\alpha} via the alfa profile log-likelihood.
}

\usage{
alfa.profile(x, a = seq(-1, 1, by = 0.01))
}

\arguments{
  \item{x}{
  A matrix with the compositional data. Zero values are not allowed.
}
  \item{a}{
  A grid of values of \eqn{\alpha}.
}
}

\details{
For every value of \eqn{\alpha} the normal likelihood (see the refernece) is computed. At the end, the plot of the values is constructed.
}

\value{
A list including:
  \item{res}{
The chosen value of \eqn{\alpha}, the corresponding log-likelihood value and the log-likelihood when \eqn{\alpha=0}.
}
  \item{ci}{
An asympotic 95\% confidence interval computed from the log-likelihood ratio test.
}
}

\references{
Tsagris M.T., Preston S. and Wood A.T.A. (2011). A data-based power transformation for compositional data. In Proceedings of the 4th Compositional Data Analysis Workshop, Girona, Spain.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{alfa.tune}, \link{alfa}, \link{alfainv}
}
}

\examples{
library(MASS)
x <- iris[, 1:4]
alfa.tune(x)
alfa.profile(x)
}

\keyword{ \eqn{\alpha}-transformation }
\keyword{ profile log-likelihood }
