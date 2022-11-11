\name{Contour plot of the alpha multivariate normal in S^2}
\alias{alfa.contour}
\title{
Contour plot of the \eqn{\alpha} multivariate normal in \eqn{S^2}
}

\description{
Contour plot of the \eqn{\alpha} multivariate normal in \eqn{S^2}.
}

\usage{
alfa.contour(m, s, a, n = 100, x = NULL, cont.line = FALSE)
}

\arguments{
\item{m}{
The mean vector of the \eqn{\alpha} multivariate normal model.
}
\item{s}{
The covariance matrix of the \eqn{\alpha} multivariate normal model.
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
\item{cont.line}{
Do you want the contour lines to appear? If yes, set this TRUE.
}
}

\details{
The \eqn{\alpha}-transformation is applied to the compositional data and then for a grid of points within the 2-dimensional simplex, the density of the \eqn{\alpha} multivariate normal is calculated and the contours are plotted.
}

\value{
The contour plot of the \eqn{\alpha} multivariate normal appears.
}

\references{
Tsagris M. and Stewart C. (2022). A Review of Flexible Transformations for Modeling Compositional Data. In Advances and Innovations in Statistics and Data Science, pp. 225--234.
https://link.springer.com/chapter/10.1007/978-3-031-08329-7_10

Tsagris M.T., Preston S. and Wood A.T.A. (2011). A data-based power transformation for compositional data.
In Proceedings of the 4th Compositional Data Analysis Workshop, Girona, Spain.
https://arxiv.org/pdf/1106.1451.pdf
}

\author{
Michail Tsagris and Christos Adam.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} and Christos Adam \email{pada4m4@gmail.com}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{folded.contour}, \link{compnorm.contour}, \link{diri.contour}, \link{mix.compnorm.contour}, \link{bivt.contour}, \link{skewnorm.contour}
}
}

\examples{
x <- as.matrix(iris[, 1:3])
x <- x / rowSums(x)
a <- a.est(x)$best
m <- colMeans(alfa(x, a)$aff)
s <- cov(alfa(x, a)$aff)
alfa.contour(m, s, a)
}

\keyword{ contour plot }
\keyword{ bivariate normal distribution }
