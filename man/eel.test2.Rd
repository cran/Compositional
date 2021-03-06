\name{Exponential empirical likelihood hypothesis testing for two mean vectors}
\alias{eel.test2}
\title{
Exponential empirical likelihood hypothesis testing for two mean vectors
}

\description{
Exponential empirical likelihood hypothesis testing for two mean vectors.
}

\usage{
eel.test2(y1, y2, tol = 1e-07, R = 0, graph = FALSE)
}

\arguments{
\item{y1}{
A matrix containing the Euclidean data of the first group.
}
\item{y2}{
A matrix containing the Euclidean data of the second group.
}
\item{tol}{
The tolerance level used to terminate the Newton-Raphson algorithm.
}
\item{R}{
If R is 0, the classical chi-square distribution is used, if R = 1, the corrected chi-square distribution (James, 1954) is used and if R = 2, 
the modified F distribution (Krishnamoorthy and Yanping, 2006) is used. If R is greater than 3 bootstrap calibration is performed.
}
\item{graph}{
A boolean variable which is taken into consideration only when bootstrap calibration is performed. IF TRUE the histogram of the bootstrap test statistic values is plotted.
}
}

\details{
Exponential empirical likelihood is a non parametric hypothesis testing procedure for one sample. The generalization to two (or more samples) is via 
searching for the mean vector that minimises the sum of the two test statistics.
}

\value{
A list including:
\item{test}{
The empirical likelihood test statistic value.
}
\item{modif.test}{
The modified test statistic, either via the chi-square or the F distribution.
}
\item{dof}{
The degrees of freedom of the chi-square or the F distribution.
}
\item{pvalue}{
The asymptotic or the bootstrap p-value.
}
\item{mu}{
The estimated common mean vector.
}
\item{runtime}{
The runtime of the bootstrap calibration.
}
}

\references{
Jing Bing-Yi and Andrew TA Wood (1996). Exponential empirical likelihood is not Bartlett correctable. Annals of Statistics 24(1): 365-369.

G.S. James (1954). Tests of Linear Hypothese in Univariate and Multivariate Analysis
when the Ratios of the Population Variances are Unknown. Biometrika, 41(1/2): 19-43

Krishnamoorthy K. and Yanping Xia (2006).  On Selecting Tests for Equality of Two Normal Mean Vectors.
Multivariate Behavioral Research 41(4): 533-548.

Owen A. B. (2001). Empirical likelihood. Chapman and Hall/CRC Press.

Amaral G.J.A., Dryden I.L. and Wood A.T.A. (2007). Pivotal bootstrap methods for k-sample problems in directional statistics and shape analysis. 
Journal of the American Statistical Association 102(478): 695-707.

Preston S.P. and Wood A.T.A. (2010). Two-Sample Bootstrap Hypothesis Tests for Three-Dimensional Labelled Landmark Data. Scandinavian Journal of Statistics 37(4): 568-587.

Tsagris M., Preston S. and Wood A.T.A. (2017). Nonparametric hypothesis testing for equality of means on the
simplex. Journal of Statistical Computation and Simulation, 87(2): 406-422.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} and Giorgos Athineou <gioathineou@gmail.com>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{el.test2}, \link{maovjames}, \link{maov}, \link{hotel2T2}, \link{james}, \link{comp.test}
}
}

\examples{
y1 = as.matrix(iris[1:25, 1:4])
y2 = as.matrix(iris[26:50, 1:4])
eel.test2(y1, y2)
eel.test2(y1, y2 )
eel.test2( y1, y2 )
}

\keyword{ Multivariate hypothesis testing }
\keyword{ non parametric test }
