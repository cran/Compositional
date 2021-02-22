\name{The Box-Cox transformation applied to ratios of components}
\alias{bc}
\title{
The Box-Cox transformation applied to ratios of components
}

\description{
The Box-Cox transformation applied to ratios of components.
}

\usage{
bc(x, lambda)
}

\arguments{
\item{x}{
A matrix with the compositional data. The first component must be zero values free.
}
\item{lambda}{
The value of the power transformation, it has to be between -1 and 1. If zero values are present it has to
be greater than 0. If \eqn{\lambda=0} the additive log-ratio transformation (\code{\link{alr}}) is applied.
}
}

\details{
The Box-Cox transformation applied to ratios of components, as described in Aitchison (1986) is applied.
}

\value{
A matrix with the transformed data.
}

\references{
Aitchison J. (1986). The statistical analysis of compositional data. Chapman & Hall.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{alr}, \link{fp}, \link{green}, \link{alfa}
}
}

\examples{
library(MASS)
x <- as.matrix(fgl[, 2:9])
x <- x / rowSums(x)
y1 <- bc(x, 0.2)
y2 <- bc(x, 0)
rbind( colMeans(y1), colMeans(y2) )
rowSums(y1)
rowSums(y2)
}

