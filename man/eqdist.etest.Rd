\name{Energy test of equality of two distributions}
\alias{eqdist.etest}
\title{
Energy test of equality of two distributions 
}

\description{
Energy test of equality of two distributions.
}

\usage{
eqdist.etest(x, y, R = 999)
}

\arguments{
\item{x}{
A matrix with the data of the first sample. 
}
\item{y}{
A matrix with the data of the second sample. 
}
\item{R}{
The number of permutations to apply in order to compute the approximate p-value.
}
}

\details{
The energy distance of equality of two distributions is applied. The main advantage of this
implementation is that it is light-weight, memory saving, however it works for two distributions only.
}

\value{
The permutation based p-value of the energy test.
}

\references{
Szekely, G. J. and Rizzo, M. L. (2004) Testing for Equal Distributions in
High Dimension. InterStat, November (5).

Szekely, G. J. (2000) Technical Report 03-05: E-statistics: Energy of
Statistical Samples. Department of Mathematics and Statistics,
Bowling Green State University.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

\seealso{
\code{\link{aeqdist.etest}, \link{acor}, \link{acor.tune}, \link{alfa}
}
}

\examples{
x <- as.matrix(iris[1:50, 1:4])
y <- as.matrix(iris[51:100, 1:4])
eqdist.etest(x, y)
}

