\name{Contour plot of the alpha-folded model in S^2}
\alias{fold.contour}
\title{
Contour plot of the alpha-folded model in \eqn{S^2}
}

\description{
Contour plot of the alpha-folded model in \eqn{S^2}.
}

\usage{
fold.contour(m, s, p, a, n = 100, x = NULL)
}

\arguments{
\item{m}{
The mean vector of the folded model.
}
\item{s}{
The covariance matrix of the folded model.
}
\item{p}{
The probability inside the simplex of the folded model.
}
\item{a}{
The value of a for the \eqn{\alpha}-transformation.
}
\item{n}{
The number of grid points to consider over which the density is calculated.
}
\item{x}{
This is either NULL (no data) or contains a 3 column matrix with compositional data.
}
}

\details{
The \eqn{alpha}-transformation is applied to the compositional data and then for a grid of points within the 2-dimensional simplex the folded model's density is calculated and the contours are plotted.
}

\value{
The contour plot of the folded model appears.
}

\references{
Tsagris Michail and Stewart Connie (2020). A folded model for compositional data analysis.
Australian and New Zealand Journal of Statistics, 62(2): 249-277.
https://arxiv.org/pdf/1802.07330.pdf
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{diri.contour}, \link{mixnorm.contour}, \link{bivt.contour}, \link{skewnorm.contour}
}
}

\examples{
x <- as.matrix(iris[, 1:3])
x <- x / rowSums(x)
a <- a.est(x)$best
mod <- alpha.mle(x, a)
fold.contour(mod$mu, mod$su, mod$p, a)
}

\keyword{ contour plot }
\keyword{ bivariate normal distribution }
