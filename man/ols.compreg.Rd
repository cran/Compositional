\name{ols.compreg}
\alias{ols.compreg}
\title{
Non linear least squares regression
}

\description{
Non linear least squares regression.
}

\usage{
ols.compreg(y, x, B = 1, ncores = 1, xnew = NULL)
}

\arguments{
  \item{y}{
  A matrix with the compositional data (dependent variable). Zero values are allowed.
}
  \item{x}{
  The predictor variable(s), they have to be continuous.
}
  \item{B}{
  If B is greater than 1 bootstrap estimates of the standard error are returned. If B=1, no standard errors are returned.
}
  \item{ncores}{
If ncores is 2 or more parallel computing is performed. This is to be used for the case of bootstrap. If B=1, this is not taken into consideration.
}
  \item{xnew}{
  If you have new data use it, otherwise leave it NULL.
}
}

\details{
The ordinary least squares between the observed and the fitted compositional data is adopted as the objective function.
This involves numerical optimisation since the relationship is non linear. There is no log-likelihood.
}

\value{
A list including:
  \item{runtime}{
  The time required by the regression.
}
  \item{beta}{
The beta coefficients.
}
  \item{seb}{
The standard error of the beta coefficients, if bootstrap is chosen, i.e. if B > 1.
}
  \item{est}{
The fitted or the predicted values (if xnew is not NULL).
}
}

\references{
Murteira, Jose MR, and Joaquim JS Ramalho. Regression analysis of multivariate fractional data. Econometric Reviews (To appear).
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{diri.reg}, \link{esov.compreg}, \link{kl.compreg}, \link{comp.reg}, \link{comp.reg}, \link{alfa.reg}
}
}

\examples{
library(MASS)
x <- fgl[, 1]
y <- fgl[, 2:9]
mod1 <- ols.compreg(y, x, B = 1, ncores = 1)
mod2 <- esov.compreg(y, x, B = 1, ncores = 1)
}

\keyword{ ordinary least squares }
\keyword{ multivariate regression }
