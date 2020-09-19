\name{The alpha-distance}
\alias{alfadist}
\alias{alfadista}
\title{
The \eqn{\alpha}-distance
}

\description{
This is the Euclidean (or Manhattan) distance after the \eqn{\alpha}-transformation has been applied.
}

\usage{
alfadist(x, a, type = "euclidean", square = FALSE)
alfadista(xnew, x, a, type = "euclidean", square = FALSE)

}

\arguments{
\item{xnew}{
A matrix or a vector with new compositional data.
}
\item{x}{
A matrix with the compositional data.
}
\item{a}{
The value of the power transformation, it has to be between -1 and 1. If zero values are present it has to be greater than 0. If \eqn{\alpha=0}, the isometric log-ratio transformation is applied.
}
\item{type}{
Which type distance do you want to calculate after the \eqn{\alpha}-transformation, "euclidean", or "manhattan".
}
\item{square}{
In the case of the Euclidean distance, you can choose to return the squared distance by setting this TRUE.
}
}

\details{
The \eqn{\alpha}-transformation is applied to the compositional data first and then the Euclidean or the Manhattan
distance is calculated.
}

\value{
For "alfadist" a matrix including the pairwise distances of all observations or the distances between xnew and x. For "alfadista" a matrix including the pairwise distances of all observations or the distances between xnew and x.
}


\references{
Tsagris M.T., Preston S. and Wood A.T.A. (2016). Improved classification for compositional data using the
\eqn{\alpha}-transformation. Journal of Classification. 33(2):243--261.
https://arxiv.org/pdf/1506.04976v2.pdf

Tsagris M.T., Preston S. and Wood A.T.A. (2011). A data-based power transformation for compositional data.
In Proceedings of the 4th Compositional Data Analysis Workshop, Girona, Spain.
https://arxiv.org/pdf/1106.1451.pdf
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{alfa}, \link{alfainv}, \link{alfa.reg}, \link{esov}
}
}

\examples{
library(MASS)
x <- as.matrix(fgl[1:20, 2:9])
x <- x / rowSums(x)
alfadist(x, 0.1)
alfadist(x, 1)
}

\keyword{ Euclidean distance }
\keyword{ Manhattan distance }

