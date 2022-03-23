\name{Density of the folded normal distribution}
\alias{dfolded}
\title{
Density of the folded model normal distribution
}

\description{
Density of the folded model normal distribution.
}

\usage{
dfolded(x, a, p, mu, su, logged = TRUE)
}

\arguments{
\item{x}{
A vector or a matrix with compositional data. No zeros are allowed.
}
\item{a}{
The value of \eqn{\alpha}.
}
\item{p}{
The probability inside the simplex of the folded model.
}
\item{mu}{
The mean vector.
}
\item{su}{
The covariance matrix.
}
\item{logged}{
A boolean variable specifying whether the logarithm of the density values to
be returned. It is set to TRUE by default.
}
}

\details{
Density values of the folded model.
}

\value{
The density value(s).
}

\references{
Tsagris M. and Stewart C. (2020). A folded model for compositional data analysis.
Australian and New Zealand Journal of Statistics, 62(2): 249-277.
https://arxiv.org/pdf/1802.07330.pdf
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{rfolded}, \link{a.est}, \link{folded.contour}
}
}

\examples{
s <- c(0.1490676523, -0.4580818209,  0.0020395316, -0.0047446076, -0.4580818209,
1.5227259250,  0.0002596411,  0.0074836251,  0.0020395316,  0.0002596411,
0.0365384838, -0.0471448849, -0.0047446076,  0.0074836251, -0.0471448849,
0.0611442781)
s <- matrix(s, ncol = 4)
m <- c(1.715, 0.914, 0.115, 0.167)
x <- rfolded(100, m, s, 0.5)
mod <- a.est(x)
den <- dfolded(x, mod$best, mod$p, mod$mu, mod$su)
}

\keyword{ random vectors simulation }
\keyword{ Gaussian mixture models }
