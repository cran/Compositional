\name{The alpha-TFLR model}
\alias{atflr}
\title{
The \eqn{\alpha}-TFLR model for compositional responses and predictors
}

\description{
The \eqn{\alpha}-TFLR model for compositional responses and predictors.
}

\usage{
atflr(y, x, a = seq(0.1, 1, by = 0.1), xnew)
}

\arguments{
\item{y}{
A matrix with the compositional data (dependent variable). Zero values are allowed.
}
\item{x}{
A matrix with the compositional predictors. Zero values are allowed.
}
\item{a}{
A vector or a single number of values of the \eqn{\alpha}-parameter. This has to be different from zero, and it can take negative values if there are no zeros in the simplicial response (y).
}
\item{xnew}{
The new data for which predictions will be made.
}
}

\details{
This is an extension of the TFLR model that includes the \eqn{\alpha}-transformation and is intended solely for prediction purposes.
}

\value{
A list with matrices containing the predicted simplicial response values, one matrix for each value of \eqn{\alpha}.
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
\code{\link{tflr}, \link{cv.atflr}, \link{scls} }
}

\examples{
library(MASS)
set.seed(1234)
y <- rdiri(214, runif(4, 1, 3))
x <- as.matrix(fgl[, 2:9])
x <- x / rowSums(x)
mod <- ascls(y, x, a = c(0.5, 1), xnew = x)
mod
}

