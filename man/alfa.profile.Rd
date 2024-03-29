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
Tsagris M.T., Preston S. and Wood A.T.A. (2011). A data-based power transformation for compositional data.
In Proceedings of the 4th Compositional Data Analysis Workshop, Girona, Spain.
https://arxiv.org/pdf/1106.1451.pdf
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} 
and Giorgos Athineou <gioathineou@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{alfa.tune}, \link{alfa}, \link{alfainv}
}
}

\examples{
x <- as.matrix(iris[, 1:4])
x <- x / rowSums(x)
alfa.tune(x)
alfa.profile(x)
}

\keyword{ profile log-likelihood }
