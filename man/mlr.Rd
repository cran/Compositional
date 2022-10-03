\name{The multiplicative log-ratio transformation and its inverse}
\alias{mlr}
\alias{mlrinv}
\title{
The multiplicative log-ratio transformation and its inverse
}

\description{
The multiplicative log-ratio transformation and its inverse.
}

\usage{
mlr(x)
mlrinv(y)
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
The multiplicative log-ratio transformation and its inverse are applied here.
This means that no zeros are allowed.
}

\value{
A matrix with the mlr transformed data (if mlr is used) or with the compositional data (if the mlrinv is used).
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
\code{\link{alr}, \link{pivot}, \link{green}, \link{alfa}
}
}

\examples{
library(MASS)
x <- as.matrix(fgl[, 2:9])
x <- x / rowSums(x)
y <- mlr(x)
x1 <- mlrinv(y)
}

\keyword{Multiplicative log-ratio-transformation}
