\name{multivreg}
\alias{multivreg}
\title{
Multivariate linear regression
}

\description{
Multivariate linear regression.
}

\usage{
multivreg(y, x, plot = TRUE, xnew = NULL)
}

\arguments{
  \item{y}{
  A matrix with the Eucldidean (continuous) data.
}
  \item{x}{
  The predictor variable(s), they have to be continuous.
}
  \item{plot}{
  Should a plot appear or not?
}
  \item{xnew}{
  If you have new data use it, otherwise leave it NULL.
}
}

\details{
 The classical multivariate linear regression model is obtained.
}

\value{
A list including:
  \item{suma}{
  A summary as produced by \code{\link{lm}}, which inlcudes the coefficients, their standard error, t-values, p-values.
}
  \item{r.squared}{
  The value of the \eqn{R^2} for each univariate regression.
}
  \item{resid.out}{
  A vector with number indicating which vectors are potential residual outliers.
}
  \item{x.leverage}{
  A vector with number indicating which vectors are potential outliers in the predicto variables space.
}
  \item{out}{
  A vector with number indicating which vectors are potential outliers in the residuals and in the predictor variables space.
}
  \item{est}{
  The fitted values if xnew is NULL, or the predicted values otherwise.
}
}

\references{
K.V. Mardia, J.T. Kent and J.M. Bibby (1979). Multivariate Analysis. Academic Press.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{diri.reg}, \link{esov.compreg}, \link{kl.compreg}, \link{ols.compreg}, \link{comp.reg}
}
}

\examples{
library(MASS)
x <- iris[, 1:2]
y<- iris[, 3:4]
multivreg(y, x, plot = TRUE)
}

\keyword{ multivariate linear regression }
