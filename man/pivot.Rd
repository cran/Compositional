\name{The pivot coordinate transformation and its inverse}
\alias{pivot}
\alias{pivotinv}
\title{
The pivot coordinate transformation and its inverse
}

\description{
The pivot coordinate transformation and its inverse.
}

\usage{
pivot(x)
pivotinv(y)
}

\arguments{
\item{x}{
A numerical matrix with the compositional data.
}
\item{y}{
A numerical matrix with data to be closed into the simplex.
}
}

\details{
The pivot coordinate transformation and its inverse are computed. This means that no zeros are allowed.
}

\value{
A matrix with the alr transformed data (if pivot is used) or with the compositional data (if the pivotinv is used).
}

\references{
Peter Filzmoser, Karel Hron and Matthias Templ (2018). Applied
Compositional Data Analysis With Worked Examples in R (pages 49 and 51). Springer.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{alfa}, \link{alfainv}, \link{alr}, \link{green}  
}
}

\examples{
library(MASS)
x <- as.matrix(fgl[, 2:9])
x <- x / rowSums(x)
y <- pivot(x)
x1 <- alrinv(y)
}
