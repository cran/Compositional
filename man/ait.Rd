\name{The alpha-IT transformation}
\alias{ait}
\title{
The \eqn{\alpha}-IT transformation
}

\description{
The \eqn{\alpha}-IT transformation.
}

\usage{
ait(x, a, h = TRUE)
}

\arguments{
\item{x}{
A matrix with the compositional data.
}
\item{a}{
The value of the power transformation, it has to be between -1 and 1. If zero values are present it has to
be greater than 0. If \eqn{\alpha=0} the isometric log-ratio transformation is applied.
}
\item{h}{
A boolean variable. If is TRUE (default value) the multiplication with the Helmert sub-matrix will take place.
When \eqn{\alpha=0} and h = FALSE, the result is the centred log-ratio transformation (Aitchison, 1986).
In general, when h = FALSE the resulting transformation maps the data onto a singualr space. The sum of the vectors is equal to 0.
Hence, from the simplex constraint the data go to another constraint.
}
}

\details{
The \eqn{\alpha}-IT transformation is applied to the compositional data.
}

\value{
A matrix with the \eqn{\alpha}-IT transformed data.
}

\references{
Clarotto L., Allard D. and Menafoglio A. (2021). A new class of \eqn{\alpha}-transformations for
the spatial analysis of Compositional Data.
https://arxiv.org/abs/2110.07967
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{aitdist}, \link{ait.knn}, \link{alfa}, \link{green}, \link{alr}
}
}

\examples{
library(MASS)
x <- as.matrix(fgl[, 2:9])
x <- x / rowSums(x)
y1 <- ait(x, 0.2)
y2 <- ait(x, 1)
rbind( colMeans(y1), colMeans(y2) )
}

