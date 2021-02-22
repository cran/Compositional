\name{The folded power transformation}
\alias{fp}
\title{
The folded power transformation
}

\description{
The folded power transformation.
}

\usage{
fp(x, lambda)
}

\arguments{
\item{x}{
A matrix with the compositional data. Zero values are allowed.
}
\item{lambda}{
The value of the power transformation, it has to be between -1 and 1. If zero values are present it has to
be greater than 0. If \eqn{\lambda=0} the additive log-ratio transformation (\code{\link{alr}}) is applied.
If zero values are present \eqn{\lambda} must be strictly positive.
}
}

\details{
The folded power transformation is applied to the compositional data.
}

\value{
A matrix with the transformed data.
}

\references{
Atkinson, A. C. (1985). Plots, transformations and regression; an introduction to graphical methods
of diagnostic regression analysis Oxford University Press.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{alr}, \link{bc}, \link{green}, \link{alfa}
}
}

\examples{
library(MASS)
x <- as.matrix(fgl[, 2:9])
x <- x / rowSums(x)
y1 <- fp(x, 0.2)
y2 <- fp(x, 0)
rbind( colMeans(y1), colMeans(y2) )
rowSums(y1)
rowSums(y2)
}

