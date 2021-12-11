\name{Total variability}
\alias{totvar}
\title{
Total variability
}

\description{
Total variability.
}

\usage{
totvar(x, a = 0)
}

\arguments{
\item{x}{
A numerical matrix with the compositional data.
}
\item{a}{
The value of the power transformation, it has to be between -1 and 1. If zero values are present it has to be greater than 0.
If \eqn{\alpha=0} the centred log-ratio transformation is used.
}
}

\details{
The \eqn{\alpha}-transformation is applied and the sum of the variances of the transformed variables is calculated. 
This is the total variability. Aitchison (1986) used the centred log-ratio transformation, but we have extended it to 
cover more geometries, via the \eqn{\alpha}-transformation.
}

\value{
The total variability of the data in a given geometry as dictated by the value of \eqn{\alpha}.
}

\references{
Aitchison J. (1986). The statistical analysis of compositional data. Chapman & Hall.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{alfa}, \ link{alfainv,} \link{alfa.profile}, \link{alfa.tune}
}
}

\examples{
x <- as.matrix(iris[, 1:4])
x <- x / rowSums(x)
totvar(x)
}


