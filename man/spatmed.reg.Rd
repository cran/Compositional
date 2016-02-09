\name{spatmed.reg}
\alias{spatmed.reg}
\title{
Spatial median regression
}
\description{
Spatial median regression with Euclidean data.
}
\usage{
spatmed.reg(y, x, xnew = NULL)
}
\arguments{
  \item{y}{
  A matrix with the compsitional data. Zero values are not allowed.
}
  \item{x}{
  The predictor variable(s), they have to be continuous.
}
  \item{xnew}{
  If you have new data use it, otherwise leave it NULL.
}
}
\details{
The objective function is the minimization of the sum of the absolute residuals. It is the multivariate generalisation of the median regression.
This function is used by \code{link{comp.reg}}.
}
\value{
A list including:
  \item{beta}{
The beta coefficients.
}
  \item{seb}{
The standard error of the beta coefficients.
}
  \item{est}{
The fitted or the predicted values (if xnew is not NULL).
}
}
\references{
Biman Chakraborty (2003) On multivariate quantile regression. Journal of Statistical Planning and Inference http://www.stat.nus.edu.sg/export/sites/dsap/research/documents/tr01_2000.pdf
}
\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
}
%\note{
%%  ~~further notes~~
%}
\seealso{
\code{\link{multivreg}, \link{comp.reg}, \link{alfa.reg}, \link{esov.compreg}, \link{diri.reg}
}
}
\examples{
library(MASS)
x <- iris[, 3:4]
y <- iris[, 1:2]
mod1 <- spatmed.reg(y, x)
mod2 <- multivreg(y, x, plot = FALSE)
}
\keyword{ spatial median regression }
\keyword{ multivariate regression }
