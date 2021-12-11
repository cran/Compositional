\name{Divergence based regression for compositional data}
\alias{kl.compreg}
\alias{js.compreg}
\alias{tv.compreg}
\alias{symkl.compreg}
\title{
Divergence based regression for compositional data
}

\description{
Regression for compositional data based on the Kullback-Leibler the Jensen-Shannon divergence and the symmetric Kullback-Leibler divergence.
}

\usage{
kl.compreg(y, x, B = 1, ncores = 1, xnew = NULL, tol = 1e-07, maxiters = 50)
js.compreg(y, x, B = 1, ncores = 1, xnew = NULL)
tv.compreg(y, x, B = 1, ncores = 1, xnew = NULL)
symkl.compreg(y, x, B = 1, ncores = 1, xnew = NULL)
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
\item{tol}{
The tolerance value to terminate the Newton-Raphson procedure.
}
\item{maxiters}{
The maximum number of Newton-Raphson iterations.
}
}

\details{
In the kl.compreg the Kullback-Leibler divergence is adopted as the objective function. In case of problematic
convergence the "multinom" function by the "nnet" package is employed. This will obviously be slower. The
js.compreg uses the Jensen-Shannon divergence and the symkl.compreg uses the symmetric Kullback-Leibler divergence.
The tv.compreg uses the Total Variation divergence. There is no actual log-likelihood for neither regression.
}

\value{
A list including:
\item{runtime}{
The time required by the regression.
}
\item{iters}{
The number of iterations required by the Newton-Raphson in the kl.compreg function.
}
\item{loglik}{
The log-likelihood. This is actually a quasi multinomial regression. This is bascially minus the half deviance, or
\eqn{- \sum_{i=1}^ny_i\log{y_i/\hat{y}_i}}.
}
\item{be}{
The beta coefficients.
}
\item{covbe}{
The covariance matrix of the beta coefficients, if bootstrap is chosen, i.e. if B > 1.
}
\item{est}{
The fitted values of xnew if xnew is not NULL.
}
}

\references{
Murteira, Jose MR, and Joaquim JS Ramalho 2016. Regression analysis of multivariate fractional data.
Econometric Reviews 35(4): 515-552.

Tsagris, Michail (2015). A novel, divergence based, regression for compositional data.
Proceedings of the 28th Panhellenic Statistics Conference, 15-18/4/2015, Athens, Greece.
https://arxiv.org/pdf/1511.07600.pdf

Endres, D. M. and Schindelin, J. E. (2003). A new metric for probability distributions.
Information Theory, IEEE Transactions on 49, 1858-1860.

Osterreicher, F. and Vajda, I. (2003). A new class of metric divergences on probability
spaces and its applicability in statistics. Annals of the Institute of Statistical
Mathematics 55, 639-653.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} 
and Giorgos Athineou <gioathineou@gmail.com>.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{diri.reg}, \link{ols.compreg}, \link{comp.reg}
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

