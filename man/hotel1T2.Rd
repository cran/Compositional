\name{Hotelling's multivariate version of the 1 sample t-test for Euclidean data}
\alias{hotel1T2}
\title{
Hotelling's multivariate version of the 1 sample t-test for Euclidean data
}

\description{
Hotelling's test for testing one Euclidean population mean vector.
}

\usage{
hotel1T2(x, M, a = 0.05, R = 999, graph = FALSE)
}

\arguments{
\item{x}{
A matrix containing Euclidean data.
}
\item{a}{
The significance level, set to 0.05 by default.
}
\item{M}{
The hypothesized mean vector.
}
\item{R}{
If R is 1 no bootstrap calibration is performed and the classical p-value via the F distribution is returned. If R is greater than 1, the bootstrap p-value is returned.
}
\item{graph}{
A boolean variable which is taken into consideration only when bootstrap calibration is performed. IF TRUE the histogram of the bootstrap test statistic values is plotted.
}
}

\details{
Multivariate hypothesis test for a one sample mean vector. This is the multivariate analogue of the one sample t-test. The p-value can be calculated either asymptotically or via bootstrap.
}

\value{
A list including:
\item{m}{
The sample mean vector.
}
\item{info}{
The test statistic, the p-value, the critical value and the degrees of freedom of the F distribution (numerator and denominator).
This is given if no bootstrap calibration is employed.
}
\item{pvalue}{
The bootstrap p-value is bootstrap is employed.
}
\item{runtime}{
The runtime of the bootstrap calibration.
}
}

\references{
K.V. Mardia, J.T. Kent and J.M. Bibby (1979). Multivariate analysis.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{eel.test1}, \link{el.test1}, \link{james}, \link{hotel2T2}, \link{maov}, \link{el.test2}, \link{comp.test}
}
}

\examples{
x <- Rfast::rmvnorm(100, numeric(10), diag( rexp(10,0.5) ) )
hotel1T2(x, numeric(10), R = 1)
hotel1T2(x, numeric(10), R = 999, graph = TRUE)
}

\keyword{ Hypothesis test }
\keyword{ mean vector }
