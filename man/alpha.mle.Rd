\name{MLE of the folded model for a given value of alpha}
\alias{alpha.mle}
\alias{a.mle}
\title{
MLE of the folded model for a given value of \eqn{\alpha}
}

\description{
MLE of the folded model for a given value of \eqn{\alpha}.
}

\usage{
alpha.mle(x, a)
a.mle(a, x)
}

\arguments{
\item{x}{
A matrix with the compositional data. No zero vaues are allowed.
}
\item{a}{
A value of \eqn{\alpha}.
}
}

\details{
This is a function for choosing or estimating the value of \eqn{\alpha} in the folded model
(Tsagris and Stewart, 2018). It is called by \code{\link{a.est}}.
}

\value{
If "alpha.mle" is called, a list including:
\item{iters}{
The nubmer of iterations the EM algorithm required.
}
\item{loglik}{
The maximimized log-likelihood of the folded model.
}
\item{p}{
The estimated probability inside the simplex of the folded model.
}
\item{mu}{
The estimated mean vector of the folded model.
}
\item{su}{
The estimated covariance matrix of the folded model.
}

If "a.mle" is called, the log-likelihood is returned only.
}

\references{
Tsagris Michail and Stewart Connie, (2018). A folded model for compositional data analysis.
https://arxiv.org/pdf/1802.07330.pdf

Tsagris M.T., Preston S. and Wood A.T.A. (2011). A data-based power transformation for compositional data.
In Proceedings of the 4th Compositional Data Analysis Workshop, Girona, Spain.
http://arxiv.org/pdf/1106.1451.pdf
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{alfa.profile}, \link{alfa}, \link{alfainv}, \link{a.est}
}
}

\examples{
x <- as.matrix(iris[, 1:4])
x <- x / rowSums(x)
mod <- alfa.tune(x)
mod
alpha.mle(x, mod[1])
}

\keyword{ maximum likelihood estimation }
