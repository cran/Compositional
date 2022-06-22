\name{Alpha generalised correlations between two compositional datasets}
\alias{acor}
\title{
\eqn{\alpha} generalised correlations between two compositional datasets
}

\description{
\eqn{\alpha} generalised correlations between two compositional datasets.
}

\usage{
acor(y, x, a, type = "dcor")
}

\arguments{
\item{y}{
A matrix with the compositional data.
}
\item{x}{
A matrix with the compositional data.
}
\item{a}{
The value of the power transformation, it has to be between -1 and 1. If zero
values are present it has to be greater than 0. If \eqn{\alpha=0} the isometric
log-ratio transformation is applied. If more than one valuesare supplied the
distance or canonical correlation are computed for all values.
}
\item{type}{
The type of correlation to compute, the distance correlation ("edist"),
the canonical correlation ("cancor") or "both".
}
}

\details{
The \eqn{\alpha}-transformation is applied to each composition and then the distance correlation
or the canonical correlation is computed. If one value of \eqn{\alpha} is supplied the type="cancor"
will return all eigenvalues. If more than one values of \eqn{\alpha} are provided then the first
eigenvalue only will be returned.
}

\value{
A vector or a matrix depending on the length of the values of \eqn{\alpha}
and the type of the correlation to be computed.
}

\references{
G.J. Szekely, M.L. Rizzo and N. K. Bakirov (2007). Measuring and Testing
Independence by Correlation of Distances. Annals of Statistics, 35(6): 2769-2794.

Tsagris M.T., Preston S. and Wood A.T.A. (2011). A data-based power
transformation for compositional data.
In Proceedings of the 4th Compositional Data Analysis Workshop, Girona, Spain.
https://arxiv.org/pdf/1106.1451.pdf
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

\seealso{
\code{\link{acor.tune}, \link{aeqdist.etest}, \link{alfa}, \link{alfa.profile}
}
}

\examples{
y <- rdiri(30, runif(3) )
x <- rdiri(30, runif(4) )
acor(y, x, a = 0.4)
}

