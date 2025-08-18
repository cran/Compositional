\name{Cross-validation for the LASSO Kullback-Leibler divergence based regression}
\alias{cv.lasso.klcompreg}
\title{
Cross-validation for the LASSO Kullback-Leibler divergence based regression
}

\description{
Cross-validation for the LASSO Kullback-Leibler divergence based regression.
}

\usage{
cv.lasso.klcompreg(y, x, alpha = 1, type = "grouped", nfolds = 10,
folds = NULL, seed = NULL, graph = FALSE)
}

\arguments{
\item{y}{
A numerical matrix with compositional data with or without zeros.
}
\item{x}{
A matrix with the predictor variables.
}
\item{alpha}{
The elastic net mixing parameter, with \eqn{0 \leq \alpha \leq 1}. The penalty is defined as a 
weighted combination of the ridge and of the Lasso regression. When \eqn{\alpha=1} LASSO is 
applied, while \eqn{\alpha=0} yields the ridge regression.
}
\item{type}{
\bold{This information is copied from the package glmnet.}. If "grouped" then a grouped lasso 
penalty is used on the multinomial coefficients for a variable. This ensures they are all in our out 
together. The default in our case is "grouped".
}
\item{nfolds}{
The number of folds for the K-fold cross-validation, set to 10 by default.
}
\item{folds}{
If you have the list with the folds supply it here. You can also leave it NULL and it will create folds.
}
\item{seed}{
You can specify your own seed number here or leave it NULL.
}
\item{graph}{
If graph is TRUE (default value) a filled contour plot will appear.
}
}

\details{
The K-fold cross-validation is performed in order to select the optimal value for \eqn{\lambda}, 
the penalty parameter in LASSO.
}

\value{
The outcome is the same as in the R package glmnet. The extra addition is that if "graph = TRUE", 
then the plot of the cross-validated object is returned. The contains the logarithm of \eqn{\lambda} 
and the deviance. The numbers on top of the figure show the number of set of coefficients for each 
component, that are not zero.
}

\references{
Alenazi, A. A. (2022). f-divergence regression models for compositional data. 
Pakistan Journal of Statistics and Operation Research, 18(4): 867--882.

Friedman, J., Hastie, T. and Tibshirani, R. (2010) Regularization Paths for Generalized Linear 
Models via Coordinate Descent. Journal of Statistical Software, Vol. 33(1), 1-22.
}

\author{
Michail Tsagris and Abdulaziz Alenazi.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} 
and Abdulaziz Alenazi \email{a.alenazi@nbu.edu.sa}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{lasso.klcompreg}, \link{lassocoef.plot}, \link{lasso.compreg}, \link{cv.lasso.compreg}, \link{kl.compreg}
}
}

\examples{
library(MASS)
y <- rdiri( 214, runif(4, 1, 3) )
x <- as.matrix( fgl[, 2:9] )
mod <- cv.lasso.klcompreg(y, x)
}

