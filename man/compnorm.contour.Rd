\name{Contour plot of the normal distribution in S^2}
\alias{compnorm.contour}
\title{
Contour plot of the normal distribution in \eqn{S^2}
}

\description{
Contour plot of the normal distribution in \eqn{S^2}.
}

\usage{
compnorm.contour(m, s, type = "alr", n = 100, x = NULL, cont.line = FALSE)
}

\arguments{
\item{m}{
The mean vector.
}
\item{s}{
The covariance matrix.
}
\item{type}{
The type of trasformation used, either the additive log-ratio ("alr"),
the isometric log-ratio ("ilr") or the pivot coordinate ("pivot") transformation.
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
The alr or the ilr transformation is applied to the compositional data at first. Then
for a grid  of points within the 2-dimensional simplex the bivariate normal density is
calculated and the contours are plotted along with the points.
}

\value{
A ternary diagram with the points (if appear = TRUE) and the bivariate normal contour lines.
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
\code{\link{diri.contour}, \link{mix.compnorm.contour}, \link{bivt.contour}, \link{skewnorm.contour}
}
}

\examples{
x <- as.matrix(iris[, 1:3])
x <- x / rowSums(x)
y <- Compositional::alr(x)
m <- colMeans(y)
s <- cov(y)
compnorm.contour(m, s)
}

\keyword{ contour plot }
\keyword{ bivariate normal distribution }
