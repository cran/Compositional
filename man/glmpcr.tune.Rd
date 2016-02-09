\name{glmpcr.tune}
\alias{glmpcr.tune}
\title{
Tuning the principal components in the GLM
}
\description{
Tuning the number of principal components in the generalised linear models.
}
\usage{
glmpcr.tune(y, x, M = 10, maxk = 10, oiko = "binomial", seed = FALSE, ncores = 2)
}
\arguments{
  \item{y}{
  A real valued vector.
}
  \item{x}{
  The predictor variables, they have to be continuous.
}
  \item{M}{
  The number of folds in the cross validation.
}
  \item{maxk}{
  The maximum number of principal components to check.
}
  \item{oiko}{
  The type of regression to be used, "binomial" if you have binary data or "poisson" if you have count data.
}
  \item{seed}{
  If seed is TRUE the folds will always be the same.
}
  \item{ncores}{
  The number of cores to use. If more than 1, parallel computing takes place.
}
}

\details{
Cross validation is performed to select the optimal number of principal components in the GLM. This is used by \code{\link{alfapcr.tune}}.
}

\value{
A list including:
  \item{mspe}{
  A vector with the mean squared error of prediction (MSPE), each value corresponds to a number of principal components.
}
  \item{k}{
  The number of principal components which minimizes the MSPE.
}
  \item{performance}{
  The bias corrected lowest value of the MSPE along with the estimated bias via the Tibshirani and Tibshirani (2009) suggestion.
}
}

\references{
Jolliffe I.T. (2002). Principal Component Analysis.

Tibshirani and Tibshirani (2009). A bias correction for the minimum error rate in cross-validation. The Annals of Applied Statistics, 3(1):822-829.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{pcr.tune}, \link{glm.pcr}, \link{alfa.pcr}, \link{alfapcr.tune}
}
}

\examples{
library(MASS)
x <- fgl[, 2:9]
y <- rpois(214, 10)
glmpcr.tune(y, x, M=10, maxk = 20, oiko = "poisson", seed = FALSE, ncores = 1)
}

\keyword{ GLM }
\keyword{ principal components }
