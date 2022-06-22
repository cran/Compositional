\name{Column-wise MLE of some univariate distributions}
\alias{colbeta.est}
\alias{collogitnorm.est}
\alias{colunitweibull.est}
\alias{colsimplex.est}

\title{
Column-wise MLE of some univariate distributions
}

\description{
Column-wise MLE of some univariate distributions.
}

\usage{
colbeta.est(x, tol = 1e-07, maxiters = 100, parallel = FALSE)
collogitnorm.est(x)
colunitweibull.est(x, tol = 1e-07, maxiters = 100, parallel = FALSE)
colsimplex.est(x, tol = 1e-07)
}

\arguments{
\item{x}{
A numerical matrix with data. Each column refers to a different vector of observations of the same distribution.
The values must by percentages, exluding 0 and 1,
}
\item{tol}{
The tolerance value to terminate the Newton-Fisher algorithm.
}
\item{maxiters}{
The maximum number of iterations to implement.
}
\item{parallel}{
Do you want to calculations to take place in parallel? The default value is FALSE
}
}

\details{
For each column, the same distribution is fitted and its parameters and log-likelihood are computed.
}

\value{
A matrix with two or three columns. The first one or the first two contain the parameter(s) of the distribution
and the second or third column the relevant log-likelihood.
}

\references{
N.L. Johnson, S. Kotz & N. Balakrishnan (1994). Continuous Univariate Distributions, Volume 1 (2nd Edition).

N.L. Johnson, S. Kotz & N. Balakrishnan (1970). Distributions in statistics: continuous univariate distributions,
Volume 2.

J. Mazucheli, A. F. B. Menezes, L. B. Fernandes, R. P. de Oliveira & M. E. Ghitany (2020).
The unit-Weibull distribution as an alternative to the Kumaraswamy distribution for the modeling of
quantiles conditional on covariates. Journal of Applied Statistics, DOI:10.1080/02664763.2019.1657813.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{censpois.mle}, \link{gammapois.mle}
}
}

\examples{
x <- matrix( rbeta(200, 3, 4), ncol = 4 )
a <- colbeta.est(x)
}
