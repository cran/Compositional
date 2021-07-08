\name{Cross-validation for the LASSO log-ratio regression with compositional response}
\alias{cv.lasso.compreg}
\title{
Cross-validation for the LASSO log-ratio regression with compositional response
}

\description{
Cross-validation for the LASSO log-ratio regression with compositional response.
}

\usage{
cv.lasso.compreg(y, x, alpha = 1, nfolds = 10,
folds = NULL, seed = FALSE, graph = FALSE)
}

\arguments{
\item{y}{
A numerical matrix with compositional data. Zero values are not allowed as the additive log-ratio transformation (\code{\link{alr}}) is applied to the compositional response prior to implementing the LASSO algortihm.
}
\item{x}{
A matrix with the predictor variables.
}
\item{alpha}{
The elastic net mixing parameter, with \eqn{0 \leq \alpha \leq 1}. The penalty is defined as a weighted
combination of the ridge and of the Lasso regression. When \eqn{\alpha=1} LASSO is applied, while
\eqn{\alpha=0} yields the ridge regression.
}
\item{nfolds}{
The number of folds for the K-fold cross validation, set to 10 by default.
}
\item{folds}{
If you have the list with the folds supply it here. You can also leave it NULL and it will create folds.
}
\item{seed}{
If seed is TRUE the results will always be the same.
}
\item{graph}{
If graph is TRUE (default value) a filled contour plot will appear.
}
}

\details{
The K-fold cross validation is performed in order to select the optimal value for \eqn{\lambda}, the
penalty parameter in LASSO.
}

\value{
The outcome is the same as in the R package glmnet. The extra addition is that if "graph = TRUE", then the
plot of the cross-validated object is returned. The contains the logarithm of \eqn{\lambda} and the mean
squared error. The numbers on top of the figure show the number of set of coefficients for each component,
that are not zero.
}

\references{
Aitchison J. (1986). The statistical analysis of compositional data. Chapman & Hall.

Friedman, J., Hastie, T. and Tibshirani, R. (2010) Regularization Paths for Generalized Linear Models via
Coordinate Descent. Journal of Statistical Software, Vol. 33(1), 1-22.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{lasso.compreg}, \link{lasso.klcompreg}, \link{lassocoef.plot}, \link{cv.lasso.klcompreg},
\link{comp.reg}
}
}

\examples{
library(MASS)
y <- rdiri( 214, runif(4, 1, 3) )
x <- as.matrix( fgl[, 2:9] )
mod <- cv.lasso.compreg(y, x)
}

