\name{Estimation of the value of alpha in the folded model}
\alias{a.est}
\title{
Estimation of the value of \eqn{\alpha} in the folded model
}

\description{
Estimation of the value of \eqn{\alpha} in the folded model.
}

\usage{
a.est(x)
}

\arguments{
\item{x}{
A matrix with the compositional data. No zero vaues are allowed.
}
}

\details{
This is a function for choosing or estimating the value of \eqn{\alpha}
in the folded model (Tsagris and Stewart, 2020).
}

\value{
A list including:
\item{runtime}{
The runtime of the algorithm.
}
\item{best}{
The estimated optimal \eqn{\alpha} of the folded model.
}
\item{loglik}{
The maximimised log-likelihood of the folded model.
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
}

\references{
Tsagris M. and Stewart C. (2022). A Review of Flexible Transformations for Modeling Compositional Data. In Advances and Innovations in Statistics and Data Science, pp. 225--234.
https://link.springer.com/chapter/10.1007/978-3-031-08329-7_10

Tsagris M. and Stewart C. (2020). A folded model for compositional data analysis.
Australian and New Zealand Journal of Statistics, 62(2): 249-277.
https://arxiv.org/pdf/1802.07330.pdf

Tsagris M.T., Preston S. and Wood A.T.A. (2011). A data-based power transformation
for compositional data. In Proceedings of the 4th Compositional Data Analysis
Workshop, Girona, Spain.
https://arxiv.org/pdf/1106.1451.pdf
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{alfa.profile}, \link{alfa}, \link{alfainv}, \link{alpha.mle}
}
}

\examples{
x <- as.matrix(iris[, 1:4])
x <- x / rowSums(x)
alfa.tune(x)
a.est(x)
}

\keyword{ maximum likelihood estimation }
