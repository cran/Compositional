\name{Contour plot of the skew skew-normal distribution in S^2}
\alias{skewnorm.contour}
\title{
Contour plot of the skew skew-normal distribution in \eqn{S^2}
}

\description{
Contour plot of the skew skew-normal distribution in \eqn{S^2}.
}

\usage{
skewnorm.contour(x, type = "alr", n = 100, appear = TRUE, cont.line = FALSE)
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
\item{cont.line}{
Do you want the contour lines to appear? If yes, set this TRUE.
}
}

\details{
The alr or the ilr transformation is applied to the compositional data at first. Then for a grid of points within the 2-dimensional
simplex the bivariate skew skew-normal density is calculated and the contours are plotted along with the points.
}

\value{
A ternary diagram with the points (if appear = TRUE) and the bivariate skew skew-normal contour lines.
}

\references{
Azzalini A. and Valle A. D. (1996). The multivariate skew-skewnormal distribution. Biometrika 83(4): 715--726.

Aitchison J. (1986). The statistical analysis of compositional data. Chapman & Hall.
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
\code{\link{diri.contour}, \link{mix.compnorm.contour}, \link{bivt.contour}, \link{compnorm.contour}
}
}

\examples{
x <- as.matrix(iris[51:100, 1:3])
x <- x / rowSums(x)
skewnorm.contour(x)
}

\keyword{ contour plot }
\keyword{ bivariate skew skewnormal distribution }
