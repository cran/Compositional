\name{The SCLS model}
\alias{scls}
\title{
Simplicial constrained linear least squares (SCLS) for compositional responses and predictors
}

\description{
Simplicial constrained linear least squares (SCLS) for compositional responses and predictors.
}

\usage{
scls(y, x, xnew = NULL, nbcores = 4)
}

\arguments{
\item{y}{
A matrix with the compositional data (dependent variable). Zero values are allowed. It may also by a big matrix of the FBM class.
}
\item{x}{
A matrix with the compositional predictors. Zero values are allowed. It may also by a big matrix of the FBM class.
}
\item{xnew}{
If you have new data use it, otherwise leave it NULL.
}
\item{nbcores}{
The number of cores to use in the case of an FBM class (big) matrix. If you do not know how many to cores to use, you may try the command \bold{nb_cores()} from the \bold{bigparallelr} package.
}
}

\details{
The function performs least squares regression where the beta coefficients are constained to be positive and sum to 1. We were inspired by the transformation-free linear regression for compositional responses and predictors of Fiksel, Zeger and Datta (2022). Our implementation now uses quadratic programming instead of the function \code{\link{optim}}, and the solution is more accurate and extremely fast.

Big matrices, of FBM class, are now accepted.
}

\value{
A list including:
\item{mse}{
The mean squared error.
}
\item{be}{
The beta coefficients.
}
\item{est}{
The fitted of xnew if xnew is not NULL.
}
}

\references{
Tsagris. M. (2025). Constrained least squares simplicial-simplicial regression.
Statistics and Computing, 35(27).

Fiksel J., Zeger S. and Datta A. (2022). A transformation-free linear regression for compositional outcomes and predictors. Biometrics, 78(3): 974--987.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
  %%  ~~further notes~~
    %}

\seealso{
\code{\link{cv.scls}, \link{tflr}, \link{scls.indeptest}, \link{scrq}  }
}

\examples{
library(MASS)
set.seed(1234)
y <- rdiri(214, runif(4, 1, 3))
x <- as.matrix(fgl[, 2:9])
x <- x / rowSums(x)
mod <- scls(y, x)
mod
}

