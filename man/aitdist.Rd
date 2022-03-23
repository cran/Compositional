\name{The alpha-IT-distance}
\alias{aitdist}
\alias{aitdista}
\title{
The \eqn{\alpha}-IT-distance
}

\description{
This is the Euclidean (or Manhattan) distance after the
\eqn{\alpha}-IT-transformation has been applied.
}

\usage{
aitdist(x, a, type = "euclidean", square = FALSE)
aitdista(xnew, x, a, type = "euclidean", square = FALSE)

}

\arguments{
\item{xnew}{
A matrix or a vector with new compositional data.
}
\item{x}{
A matrix with the compositional data.
}
\item{a}{
The value of the power transformation, it has to be between -1 and 1.
If zero values are present it has to be greater than 0. If \eqn{\alpha=0},
the isometric log-ratio transformation is applied.
}
\item{type}{
Which type distance do you want to calculate after the
\eqn{\alpha}-transformation, "euclidean", or "manhattan".
}
\item{square}{
In the case of the Euclidean distance, you can choose to return the squared
distance by setting this TRUE.
}
}

\details{
The \eqn{\alpha}-IT-transformation is applied to the compositional data first
and then the Euclidean or the Manhattan distance is calculated.
}

\value{
For "alfadist" a matrix including the pairwise distances of all observations
or the distances between xnew and x.
For "alfadista" a matrix including the pairwise distances of all observations
or the distances between xnew and x.
}


\references{
Clarotto L., Allard D. and Menafoglio A. (2021). A new class of
\eqn{\alpha}-transformations for the spatial analysis of Compositional Data.
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
\code{\link{ait}, \link{alfadist}, \link{alfa}
}
}

\examples{
library(MASS)
x <- as.matrix(fgl[1:20, 2:9])
x <- x / rowSums(x)
aitdist(x, 0.1)
aitdist(x, 1)
}

\keyword{ Euclidean distance }
\keyword{ Manhattan distance }

