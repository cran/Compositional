\name{Non linear least squares regression for compositional data}
\alias{ols.compreg}
\title{
Non linear least squares regression for compositional data
}

\description{
Non linear least squares regression for compositional data.
}

\usage{
ols.compreg(y, x, con = TRUE, B = 1, ncores = 1, xnew = NULL)
}

\arguments{
\item{y}{
A matrix with the compositional data (dependent variable). Zero values are allowed.
}
\item{x}{
A matrix or a data frame with the predictor variable(s).
}
\item{con}{
If this is TRUE (default) then the constant term is estimated, otherwise the model includes no constant term.
}
\item{B}{
If B is greater than 1 bootstrap estimates of the standard error are returned.
If B=1, no standard errors are returned.
}
\item{ncores}{
If ncores is 2 or more parallel computing is performed. This is to be used for the
case of bootstrap. If B=1, this is not taken into consideration.
}
\item{xnew}{
If you have new data use it, otherwise leave it NULL.
}
}

\details{
The ordinary least squares between the observed and the fitted compositional data
is adopted as the objective function. This involves numerical optimization since
the relationship is non linear. There is no log-likelihood.
}

\value{
A list including:
\item{runtime}{
The time required by the regression.
}
\item{beta}{
The beta coefficients.
}
\item{covbe}{
The covariance matrix of the beta coefficients, if bootstrap is chosen, i.e. if B > 1.
}
\item{est}{
The fitted of xnew if xnew is not NULL.
}
}

\references{
Murteira, Jose MR, and Joaquim JS Ramalho 2016. Regression analysis of multivariate fractional data.
Econometric Reviews 35(4): 515-552.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{diri.reg}, \link{js.compreg}, \link{kl.compreg}, \link{comp.reg}, \link{comp.reg}, \link{alfa.reg}
}
}

\examples{
library(MASS)
x <- as.vector(fgl[, 1])
y <- as.matrix(fgl[, 2:9])
y <- y / rowSums(y)
mod1 <- ols.compreg(y, x, B = 1, ncores = 1)
mod2 <- js.compreg(y, x, B = 1, ncores = 1)
}

\keyword{ ordinary least squares }
\keyword{ multivariate regression }
