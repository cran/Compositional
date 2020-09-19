\name{Greenacre's power transformation}
\alias{green}
\title{
Greenacre's power transformation
}

\description{
Greenacre's power transformation.
}

\usage{
green(x, theta)
}

\arguments{
\item{x}{
A matrix with the compositional data.
}
\item{theta}{
The value of the power transformation, it has to be between -1 and 1. If zero values are present it has to
be greater than 0. If \eqn{\theta=0} the log transformation is applied.
}
}

\details{
Greenacre's transformation is applied to the compositional data.
}

\value{
A matrix with the power transformed data.
}

\references{
Greenacre, M. (2009). Power transformations in correspondence analysis. Computational Statistics & Data Analysis, 53(8): 3107-3116.
http://www.econ.upf.edu/~michael/work/PowerCA.pdf
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{alfa}
}
}

\examples{
library(MASS)
x <- as.matrix(fgl[, 2:9])
x <- x / rowSums(x)
y1 <- green(x, 0.1)
y2 <- green(x, 0.2)
rbind( colMeans(y1), colMeans(y2) )
}

