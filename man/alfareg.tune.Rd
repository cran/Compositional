\name{alfareg.tune}
\alias{alfareg.tune}
\title{
Tuning the value of \eqn{\alpha} in the \eqn{\alpha}-regression
}

\description{
Tuning the value of \eqn{\alpha} in the \eqn{\alpha}-regression.
}

\usage{
alfareg.tune(y, x, a = seq(0.1, 1, by = 0.1), K = 10, nc = 2, graph = TRUE)
}

\arguments{
  \item{y}{
  A matrix with the compositional data. zero values are allowed.
}
  \item{x}{
  A matrix with the continuous predictor variables.
}
  \item{a}{
  The value of the power transformation, it has to be between -1 and 1. If zero values are present it has to be greater than 0. If \eqn{\alpha=0} the isometric log-ratio transformation is applied.
}
  \item{K}{
The number of folds to split the data.
}
 \item{nc}{
The number of cores to use. IF you have a multicore computer it is advisable to use more than 1. It makes the procedure faster.
}
 \item{graph}{
If grpah is TRUE a plot of the performance for each fold along the values of \eqn{\alpha} will appear.
}
}

\details{
The \eqn{\alpha}-transformation is applied to the compositional data and the numerical otpimiation is performed for the regression, unless \eqn{\alpha=0}, where the coefficients are available in closed form.
The estimated bias correction via the Tibshirani and Tibshirani (2009) criterion is applied.
}

\value{
A plot of the estimated Kullback-Leibler divergences (multiplied by 2) along the values of \eqn{\alpha} (if graph is set to TRUE).
A list including:
  \item{kl}{
Twice the Kullback-Leibler divergence of the observed from the fitted values.
}
  \item{opt}{
The optimal value of \eqn{\alpha}.
}
  \item{value}{
The minimum value of twice the Kullback-Leibler with the estimated bias added.
}
  \item{bias}{
The stimated bias.
}
}

\references{
Tsagris Michail (2015). Regression analysis with compositional data containing zero values. Chilean Journal of Statistics, 6(2): 47-57.
http://arxiv.org/pdf/1508.01913v1.pdf

Tsagris M.T., Preston S. and Wood A.T.A. (2011). A data-based power transformation for compositional data.
In Proceedings of the 4th Compositional Data Analysis Workshop, Girona, Spain.
http://arxiv.org/pdf/1106.1451.pdf

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
\code{\link{alfa}, \link{alfa.reg}
}
}

\examples{
library(MASS)
y <- fgl[1:30, 2:4]
x <- fgl[1:30, 1]
mod <- alfareg.tune(y, x, a = seq(0.2, 0.3), K = 5, nc = 1)
}

\keyword{ \eqn{\alpha}-transformation }
\keyword{ \eqn{\alpha}-regression }
