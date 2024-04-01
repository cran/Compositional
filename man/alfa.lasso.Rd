\name{LASSO with compositional predictors using the alpha-transformation}
\alias{alfa.lasso}
\title{
LASSO with compositional predictors using the \eqn{alpha}-transformation
}

\description{
LASSO with compositional predictors using the \eqn{alpha}-transformation.
}

\usage{
alfa.lasso(y, x, a = seq(-1, 1, by = 0.1), model = "gaussian", lambda = NULL,
xnew = NULL)
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
\item{xnew}{
If you have new data use it, otherwise leave it NULL.
}
}

\details{
The function uses the glmnet package to perform LASSO penalised regression. For more details see the
function in that package.
}

\value{
A list including sublists for each value of \eqn{\alpha}:
\item{mod}{
We decided to keep the same list that is returned by glmnet. So, see the function in that package
for more information.
}
\item{est}{If you supply a matrix in the "xnew" argument this will return an array of many matrices
with the fitted values, where each matrix corresponds to each value of \eqn{\lambda}.
}
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
\code{\link{alfalasso.tune}, \link{cv.lasso.klcompreg}, \link{lasso.compreg}, \link{alfa.knn.reg}
}
}

\examples{
y <- as.matrix(iris[, 1])
x <- rdiri(150, runif(20, 2, 5) )
mod <- alfa.lasso(y, x, a = c(0, 0.5, 1))
}

