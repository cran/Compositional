\name{comp.reg}
\alias{comp.reg}
\title{
Regression with compositional data
}
\description{
Regression with compositional data.
}
\usage{
comp.reg(y, x, type = "classical", xnew = NULL)
}
\arguments{
  \item{y}{
  A matrix with the compsitional data. Zero values are not allowed.
}
  \item{x}{
  The predictor variable(s), they have to be continuous.
}
  \item{type}{
  The type of regression to be used, "classical" for standard multivariate regression, "spatial" for spatial median regression, which is also robust.
}
  \item{xnew}{
  This is by default set to NULL. If you have new data whose compositional data values you want to predict, put them here.
}
}
\details{
The additive log-ratio transformation is applied and then the chosen multivariate regression is implemented. The alr is easier to explain than the ilr and that is why the latter is avoided here.
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
The fitted values if xnew is NULL, or the predicted values otherwise.
}
}
\references{
Mardia K.V., Kent J.T., and Bibby J.M. (1979). Multivariate analysis. Academic press.

Aitchison J. (1986). The statistical analysis of compositional data. Chapman & Hall.
}
\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
}
%\note{
%%  ~~further notes~~
%}
\seealso{
\code{\link{multivreg}, \link{spatmed.reg}, \link{esov.compreg}, \link{diri.reg}
}
}
\examples{
library(MASS)
y <- iris[, 1:3]
x <- iris[, 4]
mod1 <- comp.reg(y, x)
mod2 <- comp.reg(y, x, type = "spatial")
mod1
mod2
}
\keyword{ multivariate rgression }
\keyword{ spatial median regression }
