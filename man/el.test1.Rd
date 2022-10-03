\name{Empirical likelihood for a one sample mean vector hypothesis testing}
\alias{el.test1}
\title{
Empirical likelihood for a one sample mean vector hypothesis testing
}

\description{
Empirical likelihood for a one sample mean vector hypothesis testing.
}

\usage{
el.test1(x, mu, R = 1, ncores = 1, graph = FALSE)
}

\arguments{
\item{x}{
A matrix containing Euclidean data.
}
\item{mu}{
The hypothesized mean vector.
}
\item{R}{
If R is 1 no bootstrap calibration is performed and the classical p-value via the \eqn{\chi^2} distribution is returned. If R is greater than 1, the bootstrap p-value is returned.
}
\item{ncores}{
The number of cores to use, set to 1 by default.
}
\item{graph}{
A boolean variable which is taken into consideration only when bootstrap calibration is performed. IF TRUE the histogram of the bootstrap test statistic values is plotted.
}
}

\details{
Multivariate hypothesis test for a one sample mean vector. This is a non-parametric test and it works for univariate and multivariate data. 
}

\value{
A list with the outcome of the function \code{\link{el.test}} which includes
the -2 log-likelihood ratio, the observed P-value by chi-square approximation, the final value of Lagrange multiplier \eqn{\lambda},
the gradient at the maximum, the Hessian matrix, the weights on the observations (probabilities multiplied by the sample size) and the number of iteration performed. 
In addition the runtime of the procedure is reported. In the case of bootstrap, the bootstrap p-value is also returned.
}

\references{
Owen, A. (1990). Empirical likelihood ratio confidence regions. Annals of Statistics, 18, 90-120.

Owen A. B. (2001). Empirical likelihood. Chapman and Hall/CRC Press.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} 
and Giorgos Athineou <gioathineou@gmail.com>.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{eel.test1}, \link{hotel1T2}, \link{james}, \link{hotel2T2}, \link{maov}, \link{el.test2}, \link{comp.test}
}
}

\examples{
x <- Rfast::rmvnorm(100, numeric(10), diag( rexp(10, 0.5) ) )
el.test1(x, mu = numeric(10) )
eel.test1(x, mu = numeric(10) )
}

\keyword{ Hypothesis test }
\keyword{ mean vector }
