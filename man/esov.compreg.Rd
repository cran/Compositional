\name{esov.compreg}
\alias{esov.compreg}

\title{ESOV regression}

\description{
Regression based on the ESOV divergence.
}

\usage{
esov.compreg(y, x, B = 1000, ncores = 1, xnew = NULL)
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
The ESOV metric is adopted as the objective function. This involves numerical optimisation. There is no log-likelihood.
}

\value{
A list including:
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
Michail Tsagris (2015). A novel, divergence based, regression for compositional data.
Proceedings of the 28th Panhellenic Statistics Conference. http://arxiv.org/pdf/1511.07600v1.pdf

Endres, D. M. and Schindelin, J. E. (2003). A new metric for probability distributions. Information Theory, IEEE Transactions on, 49(7):1858-1860.

Osterreicher, F. and Vajda, I. (2003). A new class of metric divergences on probability spaces and its applicability in statistics. Annals of the Institute of Statistical Mathematics, 55(3):639-653.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{diri.reg}, \link{kl.compreg}, \link{ols.compreg}, \link{comp.reg}, \link{alfa.reg}
}
}

\examples{
library(MASS)
x <- fgl[, 1]
y <- fgl[, 2:9]
mod1 <- esov.compreg(y, x, B = 1, ncores = 1)
mod2 <- kl.compreg(y, x, B = 1, ncores = 1)
mod1
mod2
}

\keyword{ ESOV metric }
\keyword{ multivariate regression }
