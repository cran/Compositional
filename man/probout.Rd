\name{Estimation of the probability left outside the simplex when using the alpha-transformation}
\alias{probout}
\title{
Estimation of the probability left outside the simplex when using the alpha-transformation
}

\description{
Estimation of the probability left outside the simplex when using the alpha-transformationn.
}

\usage{
probout(mu, su, a) 
}

\arguments{
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
When applying the \eqn{\alpha}-transformation based on a multivariate normal there might be 
probability left outside the simplex as the space of this transformation is a subspace of the 
Euclidean space. The function estimates the missing probability via Monte Carlo simulation using
40 million generated vectors.
}

\value{
The estimated probability left outside the simplex.
}

\references{
Tsagris M. and Stewart C. (2020). A folded model for compositional data analysis. 
Australian and New Zealand Journal of Statistics, 62(2): 249-277.
https://arxiv.org/pdf/1802.07330.pdf

Tsagris M.T., Preston S. and Wood A.T.A. (2011). A data-based power transformation 
for compositional data. In Proceedings of the 4th Compositional Data Analysis Workshop, 
Girona, Spain.
https://arxiv.org/pdf/1106.1451.pdf
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{alfa}, \link{alpha.mle}, \link{a.est}, \link{rfolded}
}
}

\examples{
\donttest{
s <-  c(0.1490676523, -0.4580818209,  0.0020395316, -0.0047446076, -0.4580818209,
1.5227259250,  0.0002596411,  0.0074836251,  0.0020395316,  0.0002596411,
0.0365384838, -0.0471448849, -0.0047446076,  0.0074836251, -0.0471448849,
0.0611442781)
s <- matrix(s, ncol = 4)
m <- c(1.715, 0.914, 0.115, 0.167)
probout(m, s, 0.5)
}
}
