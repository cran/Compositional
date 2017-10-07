\name{James multivariate version of the t-test}
\alias{james}
\title{
James multivariate version of the t-test
}

\description{
James test for testing the equality of two population mean vectors without assuming equality of the covariance matrices.
}

\usage{
james(y1, y2, a = 0.05, R = 999, graph = FALSE)
}

\arguments{
  \item{y1}{
  A matrix containing the Euclidean data of the first group.
}
  \item{y2}{
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
Multivariate analysis of variance without assuming equality of the covariance matrices. The p-value can be calculated either asymptotically or via bootstrap. The James test (1954) or a modification proposed by Krishnamoorthy and Yanping (2006) is implemented. The James test uses a corected chi-square distribution, whereas the modified version uses an F distribution.
}

\value{
A list including:
\item{note}{
A message informing the user about the test used.
}
\item{mesoi}{
The two mean vectors.
}
\item{info}{
The test statistic, the p-value, the correction factor and the corrected critical value of the chi-square distribution if the James test has been used or,
the test statistic, the p-value, the critical value and the degrees of freedom (numerator and denominator) of the F distribution if the modified James test has been used.
}
\item{pvalue}{
The bootstrap p-value if bootstrap is employed.
}
\item{runtime}{
The runtime of the bootstrap calibration.
}
}

\references{
G.S. James (1954). Tests of Linear Hypothese in Univariate and Multivariate Analysis
when the Ratios of the Population Variances are Unknown. Biometrika, 41(1/2): 19-43

Krishnamoorthy K. and Yanping Xia.  On Selecting Tests for Equality of Two Normal Mean Vectors (2006).
Multivariate Behavioral Research 41(4): 533-548
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <gioathineou@gmail.com>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{hotel2T2}, \link{maovjames}, \link{el}, \link{comp.test}
}
}

\examples{
james( as.matrix(iris[1:25, 1:4]), as.matrix(iris[26:50, 1:4]), R = 1 )
james( as.matrix(iris[1:25, 1:4]), as.matrix(iris[26:50, 1:4]), R = 2 )
james( as.matrix(iris[1:25, 1:4]), as.matrix(iris[26:50, 1:4]) )
hotel2T2( as.matrix(iris[1:25, 1:4]), as.matrix(iris[26:50, 1:4]) )
}

\keyword{ Multivariate analysis of variance }
\keyword{ Unequality of the covariance matrices }
