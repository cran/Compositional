\name{The Frechet mean for compositional data}
\alias{frechet}
\title{
The Frechet mean for compositional data
}

\description{
Mean vector using the \eqn{\alpha}-transformation.
}

\usage{
frechet(x, a)
}

\arguments{
\item{x}{
A matrix with the compositional data.
}
\item{a}{
The value of the power transformation, it has to be between -1 and 1.
If zero values are present it has to be greater than 0. If \eqn{\alpha=0} the isometric log-ratio transformation is
applied and the closed geometric mean is calculated.
}
}

\details{
The power transformation is applied to the compositional data and the mean vector is calculated.
Then the inverse of it is calculated and the inverse of the power transformation applied to the last vector is the Frechet mean.
}

\value{
A vector with the Frechet mean for the given value of \eqn{\alpha}.
}

\references{
Tsagris M.T., Preston S. and Wood A.T.A. (2011). A data-based power transformation for compositional data.
In Proceedings of the 4th Compositional Data Analysis Workshop, Girona, Spain.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} and Giorgos Athineou <gioathineou@gmail.com>
}

%\note{
%%  ~~further notes~~
%}
\seealso{
\code{\link{alfa}, \link{alfainv}, \link{profile}
}
}

\examples{
library(MASS)
x <- as.matrix(fgl[, 2:9])
x <- x / rowSums(x)
frechet(x, 0.2)
frechet(x, 1)
}

\keyword{ Frechet mean }
