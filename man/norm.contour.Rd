\name{norm.contour}
\alias{norm.contour}
\title{
Contour plot of the normal distribution in \eqn{S^2}
}
\description{
Contour plot of the normal distribution in \eqn{S^2}.
}
\usage{
norm.contour(x, type = "alr", n = 100, appear = TRUE)
}
\arguments{
  \item{x}{
  A matrix with the compositional data. It has to be a 3 column matrix.
}
  \item{type}{
  This is either "alr" or "ilr", corresponding to the additive and the isometric log-ratio transformation respectively.
}
  \item{n}{
  The number of grid points to consider over which the density is calculated.
}
  \item{appear}{
  Should the available data appear on the ternary plot (TRUE) or not (FALSE)?
}
}
\details{
The alr or the ilr transformation is applied to the compositional data at first. Then for a grid of points within the 2-dimensional
simplex the bivariate normal density is calculated and the contours are plotted along with the points.
}
\value{
A ternary diagram with the points (if appear = TRUE) and the bivariate normal contour lines.
}
\references{

Aitchison J. (1986). The statistical analysis of compositional data. Chapman & Hall.
}
\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
}
%\note{
%%  ~~further notes~~
%}
\seealso{
\code{\link{diri.contour}, \link{mixnorm.contour}, \link{bivt.contour}, \link{skewnorm.contour}
}
}
\examples{
x <- iris[, 1:3]
x <- x / rowSums(x)
norm.contour(x)
dev.new()
norm.contour(x, type = "ilr")
}
\keyword{ contour plot }
\keyword{ bivariate normal distribution }
