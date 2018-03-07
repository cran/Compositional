\name{Hypothesis testing for two or more compositional mean vectors}
\alias{comp.test}
\title{
Hypothesis testing for two or more compositional mean vectors
}

\description{
Hypothesis testing for two or more compositional mean vectors.
}

\usage{
comp.test(x, ina, test = "james", R = 0, ncores = 1, graph = FALSE)
}

\arguments{
\item{x}{
A matrix containing compositional data.
}
\item{ina}{
A numerical or factor variable indicating the groups of the data.
}
\item{test}{
This can take the values of "james" for James' test, "hotel" for Hotelling's test, "maov" for multivariate analysis of variance assuming equality of the 
covariance matrices, "maovjames" for multivariate analysis of variance without assuming equality of the covariance matrices. "el" for empirical likelihood 
or "eel" for exponential empirical likelihood.
}

\item{R}{
This depends upon the value of the argument "test". If the test is "maov" or "maovjames", R is not taken into consideration.
If test is "hotel", then R denotes the number of bootstrap resamples. If test is "james", then R can be 1 (chi-square distribution), 
2 ( F distribution), or more for bootstrap calibration. If test is "el", then R can be 0 (chi-square), 1 (corrected chi-sqaure), 2 (F distribution) 
or more for bootstrap calibration. See the help page of each test for more information.

}
\item{ncores}{
How many to cores to use. This is taken into consideration only if test is "el" and R is more than 2.
}
\item{graph}{
A boolean variable which is taken into consideration only when bootstrap calibration is performed. IF TRUE the histogram of the bootstrap test 
statistic values is plotted. This is taken into account only when R is greater than 2.
}
}

\details{
The idea is to apply the \eqn{\alpha}-transformation, with \eqn{\alpha=1}, to the compositional data and then use a test to compare their mean vectors. 
See the help page of each test for more information. The function is visible so you can see exactly what is going on.
}

\value{
A list including:
\item{result}{
The outcome of each test.
}
}

\references{
Tsagris M., Preston S. and Wood A.T.A. (2017). Nonparametric hypothesis testing for equality of means on the
simplex. Journal of Statistical Computation and Simulation, 87(2): 406-422.

G.S. James (1954). Tests of Linear Hypothese in Univariate and Multivariate Analysis
when the Ratios of the Population Variances are Unknown. Biometrika, 41(1/2): 19-43

Krishnamoorthy K. and Yanping Xia (2006).  On Selecting Tests for Equality of Two Normal Mean Vectors.
Multivariate Behavioral Research 41(4): 533-548.

Owen A. B. (2001). Empirical likelihood. Chapman and Hall/CRC Press.

Owen A.B. (1988). Empirical likelihood ratio confidence intervals for a single functional. Biometrika 75(2): 237-249.

Amaral G.J.A., Dryden I.L. and Wood A.T.A. (2007). Pivotal bootstrap methods for k-sample problems in directional statistics and shape analysis. 
Journal of the American Statistical Association 102(478): 695-707.

Preston S.P. and Wood A.T.A. (2010). Two-Sample Bootstrap Hypothesis Tests for Three-Dimensional Labelled Landmark Data. Scandinavian Journal of Statistics 37(4): 568-587.

Jing Bing-Yi and Andrew TA Wood (1996). Exponential empirical likelihood is not Bartlett correctable. Annals of Statistics 24(1): 365-369.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <gioathineou@gmail.com>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{maovjames}, \link{maov}, \link{hotel2T2}, \link{el.test2}
}
}

\examples{
ina <- rep(1:2, each = 50)
comp.test( as.matrix(iris[1:100, 1:4]), ina, test = "james", R = 0 )
comp.test( as.matrix(iris[1:100, 1:4]), ina, test = "hotel", R = 0 )
comp.test( as.matrix(iris[1:100, 1:4]), ina, test = "el", R = 0 )
comp.test( as.matrix(iris[1:100, 1:4]), ina, test = "eel", R = 0 )
}

\keyword{ Multivariate hypothesis testing }
\keyword{ non parametric test }
\keyword{ compositional data }
