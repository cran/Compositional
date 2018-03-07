\name{Empirical likelihood hypothesis testing for two mean vectors}
\alias{el.test2}
\title{
Empirical likelihood hypothesis testing for two mean vectors
}

\description{
Empirical likelihood hypothesis testing for two mean vectors.
}

\usage{
el.test2(y1, y2, R = 0, ncores = 1, graph = FALSE)
}

\arguments{
\item{y1}{
A matrix containing the Euclidean data of the first group.
}
\item{y2}{
A matrix containing the Euclidean data of the second group.
}
\item{R}{
If R is 0, the classical chi-square distribution is used, if R = 1, the corrected chi-square distribution (James, 1954) is used and 
if R = 2, the modified F distribution (Krishnamoorthy and Yanping, 2006) is used. If R is greater than 3 bootstrap calibration is performed.
}
\item{ncores}{
How many to cores to use.
}
\item{graph}{
A boolean variable which is taken into consideration only when bootstrap calibration is performed. IF TRUE the histogram of the bootstrap test 
statistic values is plotted.
}
}

\details{
Empirical likelihood is a non parametric hypothesis testing procedure for one sample. The generalization to two (or more samples) is via searching
for the mean vector that minimizes the sum of the two test statistics.
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
Thre degrees of freedom of the chi-square or the F distribution.
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
G.S. James (1954). Tests of Linear Hypotheses in Univariate and Multivariate Analysis
when the Ratios of the Population Variances are Unknown. Biometrika, 41(1/2): 19-43

Krishnamoorthy K. and Yanping Xia (2006).  On Selecting Tests for Equality of Two Normal Mean Vectors.
Multivariate Behavioral Research 41(4): 533-548.

Owen A. B. (2001). Empirical likelihood. Chapman and Hall/CRC Press.

Owen A.B. (1988). Empirical likelihood ratio confidence intervals for a single functional. Biometrika 75(2): 237-249.

Amaral G.J.A., Dryden I.L. and Wood A.T.A. (2007). Pivotal bootstrap methods for k-sample problems in directional statistics and shape analysis. 
Journal of the American Statistical Association 102(478): 695-707.

Preston S.P. and Wood A.T.A. (2010). Two-Sample Bootstrap Hypothesis Tests for Three-Dimensional Labelled Landmark Data. Scandinavian Journal of Statistics 37(4): 568-587.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <gioathineou@gmail.com>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{eel.test2}, \link{maovjames}, \link{maov}, \link{hotel2T2}, \link{james}, \link{comp.test}
}
}

\examples{
el.test2( y1 = as.matrix(iris[1:25, 1:4]), y2 = as.matrix(iris[26:50, 1:4]), R = 0 )
el.test2( y1 = as.matrix(iris[1:25, 1:4]), y2 = as.matrix(iris[26:50, 1:4]), R = 1 )
el.test2( y1 =as.matrix(iris[1:25, 1:4]), y2 = as.matrix(iris[26:50, 1:4]), R = 2 )
}

\keyword{ Multivariate hypothesis testing }
\keyword{ Non parametric test }
