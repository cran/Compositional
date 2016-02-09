\name{fast.alfa}
\alias{fast.alfa}
\title{
Fast estimation of the value of \eqn{\alpha}
}

\description{
Fast estimation of the value of \eqn{\alpha}.
}

\usage{
fast.alfa(x, B = 1, ncores = 1)
}

\arguments{
  \item{x}{
  A matrix with the compositional data. No zero vaues are allowed.
}
  \item{B}{
  If no (bootstrap based) confidence intervals should be returned this should be 1 and more than 1 otherwise.
}
  \item{ncores}{
  If ncores is greater than 1 parallel computing is performed.
}
}

\details{
This is a faster function than \code{\link{profile}} for choosing the value of \eqn{\alpha}.
}

\value{
A vector with the best alpha, the maximised log-likelihood and the log-likelihood at \eqn{\alpha=0}, when B = 1 (no bootstrap) or
a list including:
  \item{param}{
The best alpha and the value of the log-likelihod, along with the 95\% bootstrap based confidence intervals if B is greater than 1.
}
  \item{message}{
A message with some information about the histogram if B is greater than 1.
}
}

\references{
Tsagris M.T., Preston S. and Wood A.T.A. (2011). A data-based power transformation for compositional data. In Proceedings of the 4th Compositional Data Analysis Workshop, Girona, Spain. http://arxiv.org/pdf/1106.1451.pdf
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{alfa.profile}, \link{alfa}, \link{alfainv}
}
}

\examples{
library(MASS)
x <- iris[, 1:4]
fast.alfa(x)
alfa.profile(x)
}

\keyword{ \eqn{\alpha}-transformation }
\keyword{ maximum log-likelihood estimation }
