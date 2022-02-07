\name{Random values generation from some univariate distributions defined on the (0,1) interval}
\alias{rbeta1}
\alias{runitweibull}
\alias{rlogitnorm}
\title{
Random values generation from some univariate distributions defined on the \eqn{(0,1)} interval
}

\description{
Random values generation from some univariate distributions defined on the \eqn{(0,1)} interval.
}

\usage{
rbeta1(n, a)
runitweibull(n, a, b)
rlogitnorm(n, m, s, fast = FALSE)
}

\arguments{
\item{n}{
The sample size, a numerical value.
}
\item{a}{
The shape parameter of the beta distribution. In the case of the unit Weibull, this is the shape parameter.
}
\item{b}{
This is the scale parameter for the unit Weibull distribution.
}
\item{m}{
The mean of the univariate normal in \eqn{R}.
}
\item{s}{
The standard deviation of the univariate normal in \eqn{R}.
}
\item{fast}{
If you want a faster generation set this equal to TRUE. This will use the Rnorm() function from
the Rfast package. However, the speed is only observable if you want to simulate at least 500
(this number may vary among computers) observations. The larger the sample size the higher the speed-up.
}
}

\details{
The function genrates random values from the Be(a, 1), the unit Weibull or the univariate logistic
normal distribution.
}

\value{
A vector with the simulated data.
}

%\references{
%Tyler D. E. (1987). Statistical analysis for the angular central Gaussian distribution on the sphere. Biometrika 74(3): 579-589.
%}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{beta.est}, \link{colbeta.est}, \link{rdiri}
}
}

\examples{
x <- rbeta1(100, 3)
}

