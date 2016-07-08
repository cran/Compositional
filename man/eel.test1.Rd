\name{Exponential empirical likelihood for a one sample mean vector hypothesis testing}
\alias{eel.test1}
\title{
Exponential empirical likelihood for a one sample mean vector hypothesis testing
}

\description{
Exponential empirical likelihood for a one sample mean vector hypothesis testing.
}

\usage{
eel.test1(x, mu, tol = 1e-06, R = 1)
}

\arguments{
\item{x}{
A matrix containing Euclidean data.
}
\item{mu}{
The hypothesized mean vector.
}
\item{tol}{
The tolerance value used to stop the Newton-Raphson algorithm.
}
\item{R}{
The number of bootstrap samples used to calculate the p-value. If R = 1 (default value), no bootstrap calibration is performed
}
}

\details{
Multivariate hypothesis test for a one sample mean vector. This is a non parametric test and it works for univariate and multivariate data. The p-value is currently computed only asymptotically (no bootstrap calibration at the moment).
}

\value{
A list including:
\item{p}{
The estimated probabiities.
}
\item{lambda}{
The value of the Lagrangian parameter \eqn{\lambda}.
}
\item{iter}{
The number of iterations required by the newton-Raphson algorithm.
}
\item{info}{
The value of the log-likelihood ratio test statistic along with its corresponding p-value.
}
\item{runtime}{
The runtime of the process.
}
}

\references{
Jing Bing-Yi and Andrew TA Wood (1996). Exponential empirical likelihood is not Bartlett correctable. Annals of Statistics 24(1): 365-369.

Owen A. B. (2001). Empirical likelihood. Chapman and Hall/CRC Press.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{el.test1}, \link{hotel1T2}, \link{james}, \link{hotel2T2}, \link{maov}, \link{el.test2}, \link{comp.test}
}
}

\examples{
x <- MASS::mvrnorm(100, numeric(10), diag( rexp(10, 0.5) ) )
eel.test1(x, numeric(10) )
el.test1(x, numeric(10) )
}

\keyword{ Hypothesis test }
\keyword{ Mean vector }
