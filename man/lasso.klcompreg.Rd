\name{LASSO Kullback-Leibler divergence based regression}
\alias{lasso.klcompreg}
\title{
LASSO Kullback-Leibler divergence based regression
}

\description{
LASSO Kullback-Leibler divergence based regression.
}

\usage{
lasso.klcompreg(y, x, alpha = 1, lambda = NULL,
nlambda = 100, type = "grouped", xnew = NULL)
}

\arguments{
\item{y}{
A numerical matrix with compositional data. Zero values are allowed.
}
\item{x}{
A numerical matrix containing the predictor variables.
}
\item{alpha}{
The elastic net mixing parameter, with \eqn{0 \leq \alpha \leq 1}. The penalty is defined as a 
weighted combination of the ridge and of the Lasso regression. When \eqn{\alpha=1} LASSO is applied, 
while \eqn{\alpha=0} yields the ridge regression.
}
\item{lambda}{
\bold{This information is copied from the package glmnet.} A user supplied lambda sequence. 
Typical usage is to have the program compute its own lambda sequence based on nlambda and lambda.min.ratio. 
Supplying a value of lambda overrides this. WARNING: use with care. Avoid supplying a single value for 
lambda (for predictions after CV use predict() instead). Supply instead a decreasing sequence of lambda values. 
glmnet relies on its warms starts for speed, and its often faster to fit a whole path than compute a single fit.
}
\item{nlambda}{
\bold{This information is copied from the package glmnet.} The number of \eqn{lambda} values, default is 100.
}
\item{type}{
\bold{This information is copied from the package glmnet.}. If "grouped" then a grouped lasso penalty is 
used on the multinomial coefficients for a variable. This ensures they are all in our out together.
 The default in our case is "grouped".
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
A list including:
\item{mod}{
We decided to keep the same list that is returned by glmnet. So, see the function in that package 
for more information.
}
\item{est}{ If you supply a matrix in the "xnew" argument this will return an array of many matrices 
with the fitted values, where each matrix corresponds to each value of \eqn{\lambda}.
}
}

\references{
Aitchison J. (1986). The statistical analysis of compositional data. Chapman & Hall.

Alenazi, A. A. (2022). f-divergence regression models for compositional data. 
Pakistan Journal of Statistics and Operation Research, 867--882.

Friedman, J., Hastie, T. and Tibshirani, R. (2010) Regularization Paths for Generalized Linear Models 
via Coordinate Descent. Journal of Statistical Software, Vol. 33(1), 1--22.
}

\author{
Michail Tsagris and Abdulaziz Alenazi.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} and 
Abdulaziz Alenazi \email{a.alenazi@nbu.edu.sa}.
}

%\note{
  %%  ~~further notes~~
    %}

\seealso{
\code{\link{lassocoef.plot}, \link{cv.lasso.klcompreg}, \link{kl.compreg}, \link{lasso.compreg}, \link{ols.compreg}, \link{alfa.pcr}, \link{alfa.knn.reg}
}
}

\examples{
y <- as.matrix(iris[, 1:4])
y <- y / rowSums(y)
x <- matrix( rnorm(150 * 30), ncol = 30 )
a <- lasso.klcompreg(y, x)
}

