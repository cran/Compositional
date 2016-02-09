\name{comp.den}
\alias{comp.den}

\title{
Estimating location and scatter parameters
}

\description{
Estimating location and scatter parameters in a robust and non robust way.
}

\usage{
comp.den(x, type = "alr", dist = "normal")
}

\arguments{
  \item{x}{
  A matrix containing the compositional data. No zero values are allowed.
}
  \item{type}{
  A boolean variable indicating the transformation to be used. Either "alr" or "ilr" corresponding to the additive or the isometric log-ratio transformation respectively.
}
  \item{dist}{
  Takes values "normal", "t", "skewnorm", "rob" and "spatial". They first three options correspond to the parameters of the normal, t and skew normal distribution respectively.  If it set to "rob" the MCD estimates are computed and if set to "spatial" the spatial median and spatial sign covariance matrix are computed.
}
}

\details{
This function calculates robust and non robust estimates of location and scatter.
}

\value{
A list including:
The mean vector and covariance matrix mainly. Other parameters are also returned depending on the value of the argument "dist".
}

\references{
P. J. Rousseeuw and K. van Driessen (1999) A fast algorithm for the minimum covariance determinant estimator. Technometrics 41, 212-223.

Mardia K.V., Kent J.T., and Bibby J.M. (1979). Multivariate analysis. Academic press.

Aitchison J. (1986). The statistical analysis of compositional data. Chapman & Hall.

T. Karkkaminen and S. Ayramo (2005). On computation of spatial median for robust data mining. Evolutionary and Deterministic Methods for Design, Optimization  and Control with Applications to Industrial and Societal Problems EUROGEN 2005.

A Durre, D Vogel, DE Tyler (2014). The spatial sign covariance matrix with unknown location.  Journal of Multivariate Analysis, 130: 107-117.

J. T. Kent, D. E. Tyler and Y. Vardi (1994) A curious likelihood identity for the multivariate t-distribution. Communications in Statistics-Simulation and Computation 23, 441-453.

Azzalini A. and Dalla Valle A. (1996). The multivariate skew-normal distribution. Biometrika 83(4): 715-726.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{spat.med}, \link{sscov}, \link{multivt}
}
}
\examples{
library(MASS)
x <- iris[, 1:4]
x <- x / rowSums(x)
comp.den(x)
comp.den(x, type = "alr", dist = "t")
comp.den(x, type = "alr", dist = "spatial")
}

\keyword{ location and scatter }
\keyword{ robust estimation }
