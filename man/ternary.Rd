\name{Ternary diagram}
\alias{ternary}

\title{
Ternary diagram
}

\description{
Ternary diagram.
}

\usage{
ternary(x, means = TRUE, pca = FALSE)
}

\arguments{
  \item{x}{
  A matrix with the compositional data.
}
  \item{means}{
  A boolean variable. Should the closed geometric mean and the arithmetic mean appear (TRUE) or not (FALSE)?.
}
  \item{pca}{
  Should the first PCA calculated Aitchison (1983) described appear? If yes, then this should be TRUE, or FALSE otherwise.
}
}

\details{
The first PCA is calcualte using the centred log-ratio transformation as Aitchison (1983, 1986) suggested. If the data contain zero values, the first PCA will not be plotted.
There are two ways to create a ternary graph. The one I used here, where eqch edge is equal to 1 and the one Aitchison (1986) uses. For every given point, the sum of the distances
from the edges is equal to 1. Zeros in the data appear with green circles in the triangle and you will also see NaN in the closed geometric mean.
}

\value{
The ternary plot and a matrix with the closed geometric and the simple arithmetic mean vector.
}

\references{
Aitchison, J. (1983). Principal component analysis of compositional data. Biometrika 70(1):57-65.

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
\code{\link{comp.den}, \link{alfa}, \link{diri.contour}, \link{comp.kerncontour}
}
}

\examples{
library(MASS)
x <- fgl[, 2:4]
ternary(x, means = FALSE)
x <- iris[, 1:3]
ternary(x, pca = TRUE)
}

\keyword{ ternary plot }
\keyword{ visualisation }
