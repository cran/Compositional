\name{Multivariate analysis of variance (James test)}
\alias{maovjames}
\title{
Multivariate analysis of variance (James test)
}

\description{
Multivariate analysis of variance without assuming equality of the covariance matrices.
}

\usage{
maovjames(x, ina, a = 0.05)
}

\arguments{
\item{x}{
A matrix containing Euclidean data.
}
\item{ina}{
A numerical or factor variable indicating the groups of the data.
}
\item{a}{
The significance level, set to 0.005 by default.
}
}

\details{
James (1954) also proposed an alternative to MANOVA  when the covariance matrices are not assumed equal. The test statistic for \eqn{k} samples is
\deqn{
J=\sum_{i=1}^k\left(\bar{{\bf x}}_i-\bar{{\bf X}}\right)^T{\bf W}_i\left(\bar{{\bf x}}_i-\bar{{\bf X}}\right),
}
where \eqn{\bar{{\bf x}}_i} and \eqn{n_i} are the sample mean vector and sample size of the \eqn{i}-th sample respectively and \eqn{{\bf W}_i=\left(\frac{{\bf S}_i}{n_i}\right)^{-1}}, where \eqn{{\bf S}_i} is the covariance matrix of the \eqn{i}-sample mean vector and \eqn{\bar{{\bf X}}} is the estimate of the common mean \eqn{\bar{{\bf X}}=\left(\sum_{i=1}^k{\bf W}_i\right)^{-1}\sum_{i=1}^k{\bf W}_i\bar{{\bf x}}_i}.

Normally one would compare the test statistic with a \eqn{\chi^2_{r,1-\alpha}}, where \eqn{r=p\left(k-1\right)} are the degrees of freedom with \eqn{k} denoting the number of groups and \eqn{p} the dimensionality of the data. There are \eqn{r} constraints (how many univariate means must be equal, so that the null hypothesis, that all the mean vectors are equal, holds true), that is where these degrees of freedom come from. James (1954) compared the test statistic with a corrected \eqn{\chi^2} distribution instead. Let
\eqn{A} and \eqn{B} be
\eqn{A= 1+\frac{1}{2r}\sum_{i=1}^k\frac{\left[\text{tr}\left({\bf I}_p-{\bf W}^{-1}{\bf W}_i\right)\right]^2}{n_i-1}} and \eqn{B= \frac{1}{r\left(r+2\right)}\sum_{i=1}^k\left\lbrace\frac{\text{tr}\left[\left({\bf I}_p-{\bf W}^{-1}{\bf W}_i\right)^2\right]}{n_i-1}+\frac{\left[\text{tr}\left({\bf I}_p-{\bf W}^{-1}{\bf W}_i\right)\right]^2}{2\left(n_i-1\right)}\right\rbrace}.

The corrected quantile of the \eqn{\chi^2} distribution is given as before by
\eqn{2h\left(\alpha\right)=\chi^2\left(A+B\chi^2\right)}.
}

\value{
A vector with the next 4 elements:
\item{test}{
The test statistic.
}
\item{correction}{
The value of the correction factor.
}
\item{corr.critical}{
The corrected critical value of the chi-square distribution.
}
\item{p-value}{
The p-value of the corrected test statistic.
}
}

\references{
James G.S. (1954). Tests of Linear Hypotheses in Univariate and Multivariate Analysis
when the Ratios of the Population Variances are Unknown. Biometrika, 41(1/2): 19--43.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{hotel2T2}, \link{james}
}
}

\examples{
maovjames( as.matrix(iris[,1:4]), iris[,5] )
}
