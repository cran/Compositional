\name{pcr.tune}
\alias{pcr.tune}
\title{
Tuning of the principal components regression
}
\description{
Tuning the number of principal components in the principal components regression.
}
\usage{
pcr.tune(y, x, M = 10, maxk = 50, seed = FALSE, ncores = 2)
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
  \item{seed}{
  If seed is TRUE the folds will always be the same. 
}
  \item{ncores}{
  The number of cores to use. If more than 1, parallel computing takes place.
}
}

\details{
Cross validation is performed to select the optimal number of principal components in the regression. This is used by \code{\link{alfapcr.tune}}.
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
\code{\link{glmpcr.tune}, \link{glm.pcr}, \link{alfa.pcr}, \link{alfapcr.tune}
}
}

\examples{
library(MASS) 
x <- fgl[, 2:9]
y <- fgl[, 1]
pcr.tune(y, x, M = 10, maxk = 50, seed = FALSE, ncores = 1)
}

\keyword{ principal components }
\keyword{ linear regression }
