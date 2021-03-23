\name{Plot of the LASSO coefficients}
\alias{lassocoef.plot}
\title{
Plot of the LASSO coefficients
}

\description{
Plot of the LASSO coefficients.
}

\usage{
lassocoef.plot(lasso, lambda = TRUE)
}
\arguments{
\item{lasso}{
An object where you have saved the result of the LASSO regression. See the examples for more details.
}
\item{lambda}{
If you want the x-axis to contain the logarithm of the penalty parameter \eqn{\log(\lambda)} set this to TRUE. Otherwise the x-axis will contain the \eqn{L_1}-norm of the coefficients.
}
}

\details{
This function plots the \eqn{L_2}-norm of the coefficients of each predictor variable versus the \eqn{\log(\lambda)} or the \eqn{L_1}-norm of the coefficients. This is the same plot as the one prpduced by the glmnet package with type.coef = "2norm".
}

\value{
A plot of the \eqn{L_2}-norm of the coefficients of each predictor variable (y-axis) versus the \eqn{L_1}-norm of all the coefficients (x-axis).
}

\references{
Friedman, J., Hastie, T. and Tibshirani, R. (2010) Regularization Paths for Generalized Linear Models via Coordinate Descent. Journal of Statistical Software, Vol. 33(1), 1-22.
}

\author{
Michail Tsagris and Abdulaziz Alenazi.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} and Abdulaziz Alenazi \email{a.alenazi@nbu.edu.sa}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{lasso.klcompreg}, \link{cv.lasso.klcompreg}, \link{kl.compreg}, \link{ols.compreg}, \link{alfa.pcr}, \link{alfa.knn.reg}
}
}

\examples{
y <- as.matrix(iris[, 1:4])
y <- y / rowSums(y)
x <- matrix( rnorm(150 * 30), ncol = 30 )
a <- lasso.klcompreg(y, x)
lassocoef.plot(a)
}

