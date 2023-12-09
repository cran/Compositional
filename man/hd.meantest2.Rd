\name{Two-sample test of high-dimensional means for compositional data}
\alias{hd.meantest2}
\title{
Two-sample test of high-dimensional means for compositional data
}

\description{
Two-sample test of high-dimensional means for compositional data.
}

\usage{
hd.meantest2(y1, y2, R = 1)
}

\arguments{
\item{y1}{
A matrix containing the compositional data of the first group.
}
\item{y2}{
A matrix containing the compositional data of the second group.
}
\item{R}{
If R is 1 no bootstrap calibration is performed and the asymptotic p-value is returned. If R is greater than 1,
the bootstrap p-value is returned.
}
}

\details{
A two sample for high dimensional mean vectors of compositional data is implemented. See references for more details.
}

\value{
A vector with the test statistic value and its associated (bootstrap) p-value.
}

\references{
Cao Y., Lin W. and Li H. (2018). Two-sample tests of high-dimensional means for compositional data.
Biometrika, 105(1): 115-132.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{comp.test}
}
}

\examples{
m <- runif(200, 10, 15)
x1 <- rdiri(100, m)
x2 <- rdiri(100, m)
hd.meantest2(x1, x2)
}

\keyword{ Multivariate analysis of variance }
\keyword{ Unequality of the covariance matrices }
