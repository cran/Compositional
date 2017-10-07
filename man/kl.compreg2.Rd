\name{Helper functions for the Kullback-Leibler regression}
\alias{kl.compreg2}
\alias{klcompreg.boot}
\title{
Helper functions for the Kullback-Leibler regression
}

\description{
Helper functions for the Kullback-Leibler regression.
}

\usage{
kl.compreg2(y, x, xnew = NULL, tol = 1e-07, maxiters = 50)
klcompreg.boot(y, x, der, der2, id, b1, n, p, d, tol = 1e-07, maxiters = 50)
}

\arguments{
\item{y}{
A matrix with the compositional data (dependent variable). Zero values are allowed. For the klcomobreg.boot the first column is removed.
}
\item{x}{
The predictor variable(s), they can be either continnuous or categorical or both. In the klcompreg.boot this is the design matrix, with the ones in the first column.
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
\item{der}{
An vector to put the first derivative there.
}
\item{der2}{
An empty matrix to put the second derivatives there, the Hessian matrix will be put here.
}
\item{id}{
A help vector with indices.
}
\item{b1}{
The matrix with the initial estimated coefficients.
}
\item{n}{
The sample size
}
\item{p}{
The number of columns of the design matrix.
}
\item{d}{
The dimensionality of the simplex, that is the number of columns of the compositional data minus 1.
}
}

\details{
These are help functions for the \code{\link{kl.compreg}} function. They are not to be called directly by the user.
}

\value{
For kl.compreg2 a list including:
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
For the klcompreg.boot a matrix with the estimated coefficients.
}

\references{
  Murteira, Jose MR, and Joaquim JS Ramalho 2016. Regression analysis of multivariate fractional data.
  Econometric Reviews 35(4): 515-552.
}

\author{
  Michail Tsagris

  R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <gioathineou@gmail.com>
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
