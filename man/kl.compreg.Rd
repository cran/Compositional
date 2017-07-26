\name{Divergence based regression for compositional data}
\alias{kl.compreg}
\alias{js.compreg}
\title{
Divergence based regression for compositional data
}

\description{
Regression for compositional data based on the Kullback-Leibler and the Jensen-Shannon divergence.
}

\usage{
kl.compreg(y, x, B = 1, ncores = 1, xnew = NULL)
js.compreg(y, x, B = 1, ncores = 1, xnew = NULL)
}

\arguments{
  \item{y}{
  A matrix with the compositional data (dependent variable). Zero values are allowed.
}
  \item{x}{
  The predictor variable(s), they can be either continnuous or categorical or both.
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
The Kullback-Leibler divergence is adopted as the objective function. This involves numerical optimisation. There is no log-likelihood.
}

\value{
A list including:
  \item{runtime}{
  The time required by the regression.
}
  \item{be}{
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
Murteira, Jose MR, and Joaquim JS Ramalho 2016. Regression analysis of multivariate fractional data. 
Econometric Reviews 35(4): 515-552.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{diri.reg}, \link{js.compreg}, \link{ols.compreg}, \link{comp.reg}
}
}

\examples{
library(MASS)
x <- as.vector(fgl[, 1])
y <- as.matrix(fgl[, 2:9])
y <- y / rowSums(y)
mod1<- kl.compreg(y, x, B = 1, ncores = 1)
mod2 <- js.compreg(y, x, B = 1, ncores = 1)
}

\keyword{ Kullback-Lebler divergence }
\keyword{ multivariate regression }
