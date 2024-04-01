\name{The TFLR model with multiple compositional predictors}
\alias{tflr2}
\title{
The TFLR model with multiple compositional predictors
}

\description{
The TFLR model with multiple compositional predictors
}

\usage{
tflr2(y, x, wei = FALSE, xnew = NULL)
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
The transformation-free linear regression for compositional responses and predictors is implemented.
The function to be minized is \eqn{-\sum_{i=1}^ny_i\log{y_i/(X_iB)}}. This is a self implementation of the function that can be found in the package codalm. This function allows for more than one simplicial predictors and offers the possibility of assigning weights to each simplicial predictor.
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
Fiksel J., Zeger S. and Datta A. (2022). A transformation-free linear regression for compositional outcomes and predictors. Biometrics, 78(3): 974--987.

Tsagris. M. (2024). Constrained least squares simplicial-simplicial regression.
https://arxiv.org/pdf/2403.19835.pdf
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
mod <- tflr2(y, x)
mod
}

