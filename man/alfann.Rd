\name{The k-nearest neighbours using the alpha-distance}
\alias{alfann}
\title{
The k-nearest neighbours using the \eqn{\alpha}-distance
}

\description{
The k-nearest neighbours using the \eqn{\alpha}-distance.
}

\usage{
alfann(xnew, x, a, k = 10, rann = FALSE)
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
\item{k}{
The number of nearest neighbours to search for.
}
\item{rann}{
If you have large scale datasets and want a faster k-NN search, you can use
kd-trees implemented in the R package "Rnanoflann". In this case you must set this
argument equal to TRUE. Note however, that in this case, the only available
distance is by default "euclidean".
}
}

\details{
The \eqn{\alpha}-transformation is applied to the compositional data first
and the indices of the k-nearest neighbours using the Euclidean distance
are returned.
}

\value{
A matrix including the indices of the nearest neighbours of each xnew from x.
}


\references{
MTsagris M., Alenazi A. and Stewart C. (2023).
Flexible non-parametric regression models for compositional response data with zeros.
Statistics and Computing, 33(106).

https://link.springer.com/article/10.1007/s11222-023-10277-5

Tsagris M.T., Preston S. and Wood A.T.A. (2011). A data-based power
transformation for compositional data.
In Proceedings of the 4th Compositional Data Analysis Workshop, Girona, Spain.

https://arxiv.org/pdf/1106.1451.pdf
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{alfa.knn}, \link{comp.nb}, \link{alfa.rda}, \link{alfa.nb},
link{aknn.reg}, \link{alfa}, \link{alfainv}
}
}

\examples{
library(MASS)
xnew <- as.matrix(fgl[1:20, 2:9])
xnew <- xnew / rowSums(xnew)
x <- as.matrix(fgl[-c(1:20), 2:9])
x <- x / rowSums(x)
b <- alfann(xnew, x, a = 0.1, k = 10)
}

\keyword{ Euclidean distance }
\keyword{ Manhattan distance }

