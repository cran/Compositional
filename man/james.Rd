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
If R is 1 no bootstrap calibration is performed and the classical p-value via the
F distribution is returned. If R is greater than 1, the bootstrap p-value is returned.
}
\item{graph}{
A boolean variable which is taken into consideration only when bootstrap calibration is performed.
If TRUE the histogram of the bootstrap test statistic values is plotted.
}
}

\details{
Here we show the modified version of the two-sample \eqn{T^2} test (function \code{\link{hotel2T2}}) in the case where the two covariances matrices cannot be assumed to be equal.

James (1954) proposed a test for linear hypotheses of the population means when the variances (or the covariance matrices) are not known. Its form for two \eqn{p}-dimensional samples is:
\deqn{
T^2_u=\left(\bar{{\bf X}}_1-\bar{{\bf X}}_2\right)^T\tilde{{\bf S}}^{-1}\left(\bar{{\bf X}}_1-\bar{{\bf X}}_2\right),
}
where
\eqn{\tilde{{\bf S}}=\tilde{{\bf S}_1}+\tilde{{\bf S}_2}=\frac{{\bf S}_1}{n_1}+\frac{{\bf S}_2}{n_2}
}.

James (1954) suggested that the test statistic is compared with \eqn{2h\left(\alpha\right)}, a corrected \eqn{\chi^2} distribution whose form is
\deqn{
2h\left(\alpha\right)=\chi^2\left(A+B\chi^2\right),
}
where
\eqn{A=1+\frac{1}{2p}\sum_{i=1}^2\frac{\left(tr \tilde{{\bf S}}^{-1}\tilde{{\bf S}_i}\right)^2}{n_i-1}} and
\eqn{B=\frac{1}{p\left(p+2\right)}\left[\sum_{i=1}^2\frac{tr\left(\tilde{{\bf S}}^{-1}\tilde{{\bf S}_i}\right)^2}{n_i-1}+\frac{1}{2}\sum_{i=1}^2\frac{\left(\text{tr} \tilde{{\bf S}}^{-1}\tilde{{\bf S}_i}\right)^2}{n_i-1} \right]}.

If you want to do bootstrap to get the p-value, then you must transform the data under the null hypothesis. The estimate of the common mean is given by Aitchison (1986)
\deqn{
\hat{\pmb{\mu}}_c =
\left(n_1{\bf S}_1^{-1}+n_2{\bf S}_2^{-1}\right)^{-1}\left(n_1{\bf S}_1^{-1}\bar{{\bf X}}_1+n_2{\bf S}_2^{-1}\bar{{\bf X}}_2\right)=
\left(\tilde{{\bf S}}_1^{-1}+\tilde{{\bf S}}_2^{-1}\right)^{-1}\left(\tilde{{\bf S}}_1^{-1}\bar{{\bf X}}_1+\tilde{{\bf S}}_2^{-1}\bar{{\bf X}}_2\right).
}

The modified Nel and van der Merwe (1986) test is based on the same quadratic form as that of James (1954) but the distribution used to compare the value of the test statistic is different.
It is shown in Krishnamoorthy and Yanping (2006) that \eqn{T^2_u \sim \frac{\nu p}{\nu-p+1}F_{p,\nu-p+1}} approximately, where
\eqn{
\nu=\frac{p+p^2}{\frac{1}{n_1}\left\lbrace \text{tr}\left[ \left( {\bf S}_1\tilde{{\bf S}} \right)^2\right]+
\text{tr}\left[ \left( {\bf S}_1\tilde{{\bf S}} \right)\right]^2 \right\rbrace +
\frac{1}{n_2}\left\lbrace \text{tr}\left[ \left( {\bf S}_2\tilde{{\bf S}}\right)^2\right]+
\text{tr}\left[ \left( {\bf S}_2\tilde{{\bf S}} \right)\right]^2 \right\rbrace }.
}

The algorithm is taken by Krishnamoorthy and Yu (2004).
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
The test statistic, the p-value, the correction factor and the corrected critical value
of the chi-square distribution if the James test has been used or, the test statistic,
the p-value, the critical value and the degrees of freedom (numerator and denominator)
of the F distribution if the modified James test has been used.
}
\item{pvalue}{
The bootstrap p-value if bootstrap is employed.
}
\item{runtime}{
The runtime of the bootstrap calibration.
}
}

\references{
Aitchison J. (1986). The statistical analysis of compositional data. Chapman & Hall.

James G.S. (1954). Tests of Linear Hypothese in Univariate and Multivariate Analysis
when the Ratios of the Population Variances are Unknown. Biometrika, 41(1/2): 19--43.

Krishnamoorthy K. and Yu J. (2004). Modified Nel and Van der Merwe test for the multivariate
Behrens-Fisher problem. Statistics & Probability Letters, 66(2): 161--169.

Krishnamoorthy K. and Yanping Xia (2006). On Selecting Tests for Equality of Two Normal Mean Vectors.
Multivariate Behavioral Research, 41(4): 533--548.

Tsagris M., Preston S. and Wood A.T.A. (2017). Nonparametric hypothesis testing for
equality of means on the simplex. Journal of Statistical Computation and Simulation, 87(2): 406--422.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{hotel2T2}, \link{maovjames}, \link{el.test2}, \link{eel.test2}
}
}

\examples{
james( as.matrix(iris[1:25, 1:4]), as.matrix(iris[26:50, 1:4]), R = 1 )
james( as.matrix(iris[1:25, 1:4]), as.matrix(iris[26:50, 1:4]), R = 2 )
james( as.matrix(iris[1:25, 1:4]), as.matrix(iris[26:50, 1:4]) )
}

