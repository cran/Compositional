\name{The additive log-ratio transformation and its inverse}
\alias{alr}
\alias{alrinv}
\title{
The additive log-ratio transformation and its inverse
}

\description{
The additive log-ratio transformation and its inverse.
}

\usage{
alr(x)
alrinv(y)
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
The additive log-ratio transformation with the first component being the commn divisor is applied. The inverse of this trnasformation is also available.
}

\value{
A matrix with the alr transformed data (if alr is used) or with the compositional data (if the alrinv is used).
}

\references{
Aitchison J. (1986). The statistical analysis of compositional data. Chapman \& Hall.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{alfa}, \ link{alfainv,} \link{alfa.profile}, \link{alfa.tune}
}
}

\examples{
library(MASS)
x <- as.matrix(fgl[, 2:9])
x <- x / rowSums(x)
y <- alr(x)
x1 <- alrinv(y)
}

\keyword{ Additive log-ratio-transformation }
