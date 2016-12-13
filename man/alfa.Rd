\name{The alpha-transformation}
\alias{alfa}
\title{
The \eqn{\alpha}-transformation
}

\description{
The \eqn{\alpha}-transformation.
}

\usage{
alfa(x, a, h = TRUE)
}

\arguments{
  \item{x}{
  A matrix with the compositional data.
}
  \item{a}{
  The value of the power transformation, it has to be between -1 and 1. If zero values are present it has to be greater than 0.
  If \eqn{\alpha=0} the isometric log-ratio transformation is applied.
}
  \item{h}{
  A boolean variable. If is TRUE (default value) the multiplication with the Helmert sub-matrix will take place.
  When \eqn{\alpha=0} and h = FALSE, the result is the centred log-ratio transformation (Aitchison, 1986).
  In general, when h = FALSE the resulting transformation maps the data onto a singualr space. The sum of the vectors is equal to 0.
  Hence, from the simplex constraint the data go to another constraint.
}
}

\details{
The \eqn{\alpha}-transformation is applied to the compositional data.
}

\value{
A list including:
  \item{sa}{
The logarithm of the Jacobian determinant of the \eqn{\alpha}-transformation. This is used in the "profile" function to speed up the computations.
}
  \item{aff}{
The \eqn{\alpha}-transformed data.
}
}

\references{
Tsagris M.T., Preston S. and Wood A.T.A. (2011). A data-based power transformation for compositional data.
In Proceedings of the 4th Compositional Data Analysis Workshop, Girona, Spain. http://arxiv.org/pdf/1106.1451.pdf

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
\code{\link{alfainv}, \link{alfa.profile}, \link{alfa.tune}
}
}

\examples{
library(MASS)
x <- as.matrix(fgl[, 2:9])
x <- x / rowSums(x)
y1 <- alfa(x, 0.2)$aff
y2 <- alfa(x, 1)$aff
rbind( colMeans(y1), colMeans(y2) )
y3 <- alfa(x, 0.2)$aff
dim(y1)  ;  dim(y3)
rowSums(y1)
rowSums(y3)
}

\keyword{ \eqn{\alpha}-transformation }
