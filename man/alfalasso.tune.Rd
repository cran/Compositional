\name{Cross-validation for LASSO with compositional predictors using the alpha-transformation}
\alias{alfalasso.tune}
\title{
Cross-validation for LASSO with compositional predictors using the \eqn{alpha}-transformation
}

\description{
Cross-validation for LASSO with compositional predictors using the \eqn{alpha}-transformation.
}

\usage{
alfalasso.tune(y, x, a = seq(-1, 1, by = 0.1), model = "gaussian", lambda = NULL,
type.measure = "mse", nfolds = 10, folds = NULL, stratified = FALSE)
}

\arguments{
\item{y}{
A numerical vector or a matrix for multinomial logistic regression.
}
\item{x}{
A numerical matrix containing the predictor variables, compositional data, where
zero values are allowed..
}
\item{a}{
A vector with a grid of values of the power transformation, it has to be between -1 and 1.
If zero values are present it has to be greater than 0. If \eqn{\alpha=0} the isometric
log-ratio transformation is applied.
}
\item{model}{
The type of the regression model, "gaussian", "binomial", "poisson", "multinomial", or "mgaussian".
}
\item{lambda}{
\bold{This information is copied from the package glmnet.} A user supplied lambda sequence.
Typical usage is to have the program compute its own lambda sequence based on nlambda and lambda.min.ratio.
Supplying a value of lambda overrides this. WARNING: use with care. Avoid supplying a single value for
lambda (for predictions after CV use predict() instead). Supply instead a decreasing sequence of lambda values.
glmnet relies on its warms starts for speed, and its often faster to fit a whole path than compute a single fit.
}
\item{type.measure}{
\bold{This information is taken from the package glmnet.} The loss function to use for cross-validation. For gaussian models this can be "mse", "deviance" for logistic and poisson regression, "class" applies to binomial and multinomial logistic regression only, and gives misclassification error. "auc" is for two-class logistic regression only, and gives The area under the ROC curve. "mse" or "mae" (mean absolute error) can be used by all models.
}
\item{nfolds}{
The number of folds. Set to 10 by default.
}
\item{folds}{
If you have the list with the folds supply it here. You can also leave it NULL
and it will create folds.
}
\item{stratified}{
Do you want the folds to be created in a stratified way? TRUE or FALSE.
}
}

\details{
The function uses the glmnet package to perform LASSO penalised regression. For more details see the
function in that package.
}

\value{
A matrix with two columns and number of rows equal to the number of \eqn{\alpha} values used. Each row contains, the optimal value of the \eqn{\lambda} penalty parameter for the LASSO and the optimal value of the loss function, for each value of \eqn{\alpha}.
}

\references{
Aitchison J. (1986). The statistical analysis of compositional data. Chapman & Hall.

Friedman, J., Hastie, T. and Tibshirani, R. (2010) Regularization Paths for Generalized Linear Models
via Coordinate Descent. Journal of Statistical Software, Vol. 33(1), 1--22.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
  %%  ~~further notes~~
    %}

\seealso{
\code{\link{alfa.lasso}, \link{cv.lasso.klcompreg}, \link{lasso.compreg}, \link{alfa.knn.reg}
}
}

\examples{
y <- iris[, 1]
x <- rdiri(150, runif(20, 2, 5) )
mod <- alfalasso.tune( y, x, a = c(0.2, 0.5, 1) )
}

