\name{Tuning the principal components with GLMs}
\alias{glmpcr.tune}
\title{
Tuning the principal components with GLMs
}

\description{
Tuning the number of principal components in the generalised linear models.
}

\usage{
glmpcr.tune(y, x, M = 10, maxk = 10,
mat = NULL, ncores = 1, graph = TRUE)
}

\arguments{
  \item{y}{
  A real valued vector.
}
  \item{x}{
  A matrix with the predictor variables, they have to be continuous.
}
  \item{M}{
  The number of folds in the cross validation.
}
  \item{maxk}{
  The maximum number of principal components to check.
}
  \item{mat}{
You can specify your own folds by giving a mat, where each column is a fold. Each column contains indices of the observations.  You can also leave it NULL and it will create folds.
}
  \item{ncores}{
  The number of cores to use. If more than 1, parallel computing will take place. It is advisable to use it if you have many observations and or many variables, otherwise it will slow down th process.
}
 \item{graph}{
If graph is TRUE a plot of the performance for each fold along the values of \eqn{\alpha} will appear.
}
}

\details{
Cross validation is performed to select the optimal number of principal components in the GLM. This is used by \code{\link{alfapcr.tune}}.
}

\value{
If graph is TRUE a plot of the performance versus the number of principal components will appear.
A list including:
  \item{msp}{
  A matrix with the mean squared error of prediction (MSPE) for every fold.
}
  \item{mpd}{
  A vector with the mean squared error of prediction (MSPE), each value corresponds to a number of principal components.
}
  \item{k}{
  The number of principal components which minimizes the MSPE.
}
  \item{performance}{
  The bias corrected lowest value of the MSPE along with the estimated bias via the Tibshirani and Tibshirani (2009) suggestion.
}
  \item{runtime}{
  The time required by the cross-validation procedure.
  }
}

\references{
Aguilera A.M., Escabias M. and Valderrama M.J. (2006). Using principal components for estimating logistic regression with high-dimensional multicollinear data. Computational Statistics & Data Analysis 50(8): 1905-1924.

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
x <- as.matrix(fgl[, 2:9])
y <- rpois(214, 10)
glmpcr.tune(y, x, M=10, maxk = 20, mat = NULL, ncores = 1)
}

\keyword{ GLM }
\keyword{ principal components }
