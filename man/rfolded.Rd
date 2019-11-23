\name{Simulation of compositional data from the folded model}
\alias{rfolded}
\title{
Simulation of compositional data from the folded model
}

\description{
Simulation of compositional data from the folded model.
}

\usage{
rfolded(n, mu, su, a)
}

\arguments{
\item{n}{
The sample size
}
\item{mu}{
The mean vector.
}
\item{su}{
The covariance matrix.
}
\item{a}{
The value of \eqn{\alpha}.
}
}

\details{
A sample from the folded model is generated.
}

\value{
A matrix with compositional data.
}

\references{
Tsagris Michail and Stewart Connie, (2018). A folded model for compositional data analysis.
https://arxiv.org/pdf/1802.07330.pdf
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{alfa}, \link{alpha.mle}, \link{a.est}
}
}

\examples{
s <-  c(0.1490676523, -0.4580818209,  0.0020395316, -0.0047446076, -0.4580818209,
1.5227259250,  0.0002596411,  0.0074836251,  0.0020395316,  0.0002596411,
0.0365384838, -0.0471448849, -0.0047446076,  0.0074836251, -0.0471448849,
0.0611442781)
s <- matrix(s, ncol = 4)
m = c(1.715, 0.914, 0.115, 0.167)
x = rfolded(100, m, s, 0.5)
a.est(x)
}

\keyword{ random vectors simulation }
\keyword{ Gaussian mixture models }
