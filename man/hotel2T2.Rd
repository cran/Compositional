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
The fist case scenario is when we assume equality of the two covariance matrices. This is called the two-sample Hotelling's \eqn{T^2} test (Mardia, Kent and Bibby, 1979, pg. 131-140) and Everitt (2005, pg. 139). The test statistic is defined as
\deqn{
T^2=\frac{n_1n_2}{n_1+n_2}\left(\bar{{\bf X}}_1- \bar{{\bf X}}_2\right)^T{\bf S}^{-1}\left(\bar{{\bf X}}_1- \bar{{\bf X}}_2\right),
}
where \eqn{\bf S} is the pooled covariance matrix calculated under the assumption of equal covariance matrices
\eqn{{\bf S}=\frac{\left(n_1-1\right){\bf S}_1+\left(n_2-1\right){\bf S}_2}{n_1+n_2-2}.}
Under \eqn{H_0} the statistic \eqn{F} given by
\deqn{
F=\frac{\left( n_1+n_2-p-1 \right)T^2}{\left(n_1+n_2-2 \right)p}
}
follows the \eqn{F} distribution with \eqn{p} and \eqn{n_1+n_2-p-1} degrees of freedom. Similar to the one-sample test, an extra argument (R) indicates whether bootstrap calibration should be used or not. If R=1, then the asymptotic theory applies, if R>1, then the bootstrap p-value will be applied and the number of re-samples is equal to R. The estimate of the common mean used in the bootstrap to transform the data under the null hypothesis the mean vector of the combined sample, of all the observations.

The built-in command \code{\link{manova}} does the same thing exactly. Try it, the asymptotic \eqn{F} test is what you have to see. In addition, this command allows for more mean vector hypothesis testing for more than two groups. I noticed this command after I had written my function and nevertheless as I mention in the introduction this document has an educational character as well.
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
Everitt B. (2005). An R and S-Plus Companion to Multivariate Analysis. Springer.

Mardia K.V., Kent J.T. and Bibby J.M. (1979). Multivariate Analysis. London: Academic
Press.

Tsagris M., Preston S. and Wood A.T.A. (2017). Nonparametric hypothesis
testing for equality of means on the simplex.
Journal of Statistical Computation and Simulation, 87(2): 406--422.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{james}, \link{el.test2}, \link{eel.test2}
}
}

\examples{
hotel2T2( as.matrix(iris[1:25, 1:4]), as.matrix(iris[26:50, 1:4]) )
hotel2T2( as.matrix(iris[1:25, 1:4]), as.matrix(iris[26:50, 1:4]), R = 1 )
}

