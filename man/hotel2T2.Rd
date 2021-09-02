\name{Hotelling's multivariate version of the 2 sample t-test for Euclidean data}
\alias{hotel2T2}
\title{
Hotelling's multivariate version of the 2 sample t-test for Euclidean data
}

\description{
Hotelling's test for testing the equality of two Euclidean population mean vectors.
}

\usage{
hotel2T2(x1, x2, a = 0.05, R = 999, graph = FALSE)
}

\arguments{
\item{x1}{
A matrix containing the Euclidean data of the first group.
}
\item{x2}{
A matrix containing the Euclidean data of the second group.
}
\item{a}{
The significance level, set to 0.05 by default.
}
\item{R}{
If R is 1 no bootstrap calibration is performed and the classical p-value via the F distribution is returned. If R is greater than 1, the bootstrap p-value is returned.
}
\item{graph}{
A boolean variable which is taken into consideration only when bootstrap calibration is performed. IF TRUE the histogram of the bootstrap test statistic values is plotted.
}
}

\details{
Multivariate analysis of variance assuming equality of the covariance matrices. The p-value can be calculated either asymptotically or via bootstrap.
}

\value{
A list including:
\item{mesoi}{
The two mean vectors.
}
\item{info}{
The test statistic, the p-value, the critical value and the degrees of freedom of the F distribution (numerator and denominator).
This is given if no bootstrap calibration is employed.
}
\item{pvalue}{
The bootstrap p-value is bootstrap is employed.
}
\item{note}{
A message informing the user that bootstrap calibration has been employed.
}
\item{runtime}{
The runtime of the bootstrap calibration.
}
}

\references{
Everitt B. (2005). An R and S-Plus Companion to Multivariate Analysis p. 139-140. Springer.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{james}, \link{maov}, \link{el.test2}, \link{comp.test}
}
}

\examples{
hotel2T2( as.matrix(iris[1:25, 1:4]), as.matrix(iris[26:50, 1:4]) )
hotel2T2( as.matrix(iris[1:25, 1:4]), as.matrix(iris[26:50, 1:4]), R = 1 )
james( as.matrix(iris[1:25, 1:4]), as.matrix(iris[26:50, 1:4]) )
james( as.matrix(iris[1:25, 1:4]), as.matrix(iris[26:50, 1:4]), R = 1 )
}

\keyword{ Multivariate analysis of variance }
\keyword{ Equality of the covariance matrices }
