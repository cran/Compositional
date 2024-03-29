\name{Contour plot of the alpha-folded model in S^2}
\alias{folded.contour}
\title{
Contour plot of the \eqn{\alpha}-folded model in \eqn{S^2}
}

\description{
Contour plot of the \eqn{\alpha}-folded model in \eqn{S^2}.
}

\usage{
folded.contour(mu, su, p, a, n = 100, x = NULL, cont.line = FALSE)
}

\arguments{
\item{mu}{
The mean vector of the folded model.
}
\item{su}{
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
\item{cont.line}{
Do you want the contour lines to appear? If yes, set this TRUE.
}
}

\details{
The \eqn{\alpha}-transformation is applied to the compositional data and then for a grid
of points within the 2-dimensional simplex the folded model's density is calculated and
the contours are plotted.
}

\value{
The contour plot of the folded model appears.
}

\references{
Tsagris M. and Stewart C. (2022). A Review of Flexible Transformations for Modeling Compositional Data. In Advances and Innovations in Statistics and Data Science, pp. 225--234.
https://link.springer.com/chapter/10.1007/978-3-031-08329-7_10

Tsagris M. and Stewart C. (2020). A folded model for compositional data analysis.
Australian and New Zealand Journal of Statistics, 62(2): 249-277.
https://arxiv.org/pdf/1802.07330.pdf
}

\author{
Michail Tsagris and Christos Adam.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} and
Christos Adam \email{pada4m4@gmail.com}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{alfa.contour}, \link{compnorm.contour}, \link{diri.contour}, \link{mix.compnorm.contour},
\link{bivt.contour}, \link{skewnorm.contour}
}
}

\examples{
x <- as.matrix(iris[, 1:3])
x <- x / rowSums(x)
a <- a.est(x)$best
mod <- alpha.mle(x, a)
folded.contour(mod$mu, mod$su, mod$p, a)
}

\keyword{ contour plot }
\keyword{ bivariate normal distribution }
