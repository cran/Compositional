\name{The SCLS model with multiple compositional predictors}
\alias{scls2}
\title{
The SCLS model with multiple compositional predictors
}

\description{
The SCLS model with multiple compositional predictors.
}

\usage{
scls2(y, x, wei = FALSE, xnew = NULL)
}

\arguments{
\item{y}{
A matrix with the compositional data (dependent variable). Zero values are allowed.
}
\item{x}{
A list of matrices with the compositional predictors. Zero values are allowed.
}
\item{wei}{
Do you want weights among the different simplicial predictors? The default is FALSE.
}
\item{xnew}{
If you have new data use it, otherwise leave it NULL.
}
}

\details{
The function performs least squares regression where the beta coefficients are constained to be positive and sum to 1. We were inspired by the transformation-free linear regression for compositional responses and predictors of Fiksel, Zeger and Datta (2020). Our implementation now uses quadratic programming instead of the function \code{\link{optim}}, and the solution is more accurate and extremely fast. This function allows for more than one simplicial predictors and offers the possibility of assigning weights to each simplicial predictor.
}

\value{
A list including:
\item{ini.mse}{
The mean squared error when all simplicial predictors carry equal weight.
}
\item{ini.be}{
The beta coefficients when all simplicial predictors carry equal weight.
}
\item{mse}{
The mean squared error when the simplicial predictors carry unequal weights.
}
\item{weights}{
The weights in a vector form. A vector of length equal to the number of rows of the matrix of coefficients.
}
\item{am}{
The vector of weights, one for each simplicia predictor. The length of the vector is equal to the number of simplicial predictors.
}
\item{est}{
The fitted of xnew if xnew is not NULL.
}
}

\references{
Tsagris. M. (2025). Constrained least squares simplicial-simplicial regression.
Statistics and Computing, 35(27).
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
  %%  ~~further notes~~
    %}

\seealso{
\code{\link{cv.scls}, \link{tflr}, \link{scls.indeptest} }
}

\examples{
library(MASS)
set.seed(1234)
y <- rdiri(214, runif(4, 1, 3))
x1 <- as.matrix(fgl[, 2:9])
x <- list()
x[[ 1 ]] <- x1 / rowSums(x1)
x[[ 2 ]] <- Compositional::rdiri(214, runif(4))
mod <- scls2(y, x)
mod
}

