\name{Energy test of equality of distributions using the alpha-transformation}
\alias{aeqdist.etest}
\title{
Energy test of equality of distributions using the \eqn{\alpha}-transformation
}

\description{
Energy test of equality of distributions using the \eqn{\alpha}-transformation.
}

\usage{
aeqdist.etest(x, sizes, a = 1, R = 999)
}

\arguments{
\item{x}{
A matrix with the compositional data with all groups stacked one under the other.
}
\item{sizes}{
A numeric vector matrix with the sample sizes.
}
\item{a}{
The value of the power transformation, it has to be between -1 and 1. If zero values are present it has to be greater than 0. If \eqn{\alpha=0} the isometric log-ratio transformation is applied. If more than one values are supplied the energy distance of equality of distributions is applied for each value of \eqn{\alpha}.
}
\item{R}{
The number of permutations to apply in order to compute the approximate p-value.
}
}

\details{
The \eqn{\alpha}-transformation is applied to each composition and then the energy distance of equality of distributions is applied for each value of \eqn{\alpha} or for the single value of \eqn{\alpha}.
}

\value{
A numerical value or a numerical vector, depending on the length of the values of \eqn{\alpha}, with the approximate p-value(s) of the energy test.
}

\references{
Szekely, G. J. and Rizzo, M. L. (2004) Testing for Equal Distributions in High Dimension,
InterStat, November (5).

Szekely, G. J. (2000) Technical Report 03-05: E-statistics: Energy of Statistical Samples,
Department of Mathematics and Statistics, Bowling Green State University.

Tsagris M.T., Preston S. and Wood A.T.A. (2011). A data-based power transformation for compositional data.
In Proceedings of the 4th Compositional Data Analysis Workshop, Girona, Spain.
https://arxiv.org/pdf/1106.1451.pdf
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

\seealso{
\code{\link{acor}, \link{acor.tune}, \link{alfa}, \link{alfa.profile}
}
}

\examples{
y <- rdiri(50, c(3, 4, 5) )
x <- rdiri(60, c(3, 4, 5) )
aeqdist.etest( rbind(x, y), c(dim(x)[1], dim(y)[1]), a = c(-1, 0, 1) )
}

