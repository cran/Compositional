\name{The transformation-free linear regression (TFLR) for compositional responses and predictors}
\alias{tflr}
\title{
Transformation-free linear regression (TFLR) for compositional responses and predictors
}

\description{
Transformation-free linear regression (TFLR) for compositional responses and predictors.
}

\usage{
tflr(y, x, xnew = NULL)
}

\arguments{
\item{y}{
A matrix with the compositional response. Zero values are allowed.
}
\item{x}{
A matrix with the compositional predictors. Zero values are in general allowed, but there can be cases when these are problematic.
}
\item{xnew}{
If you have new data use it, otherwise leave it NULL.
}
}

\details{
The transformation-free linear regression for compositional responses and predictors is implemented.
The function to be minized is \eqn{-\sum_{i=1}^ny_i\log{y_i/(X_iB)}}. This is an efficient self implementation.
}

\value{
A list including:
\item{kl}{
The Kullback-Leibler divergence between the observed and the fitted response compositional data.
}
\item{be}{
The beta coefficients.
}
\item{est}{
The fitted values of xnew if xnew is not NULL.
}
}

\references{
Fiksel J., Zeger S. and Datta A. (2022). A transformation-free linear regression for compositional outcomes and predictors. Biometrics, 78(3): 974--987.

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
\code{\link{cv.tflr}, \link{scls} \link{kl.alfapcr}
}
}

\examples{
library(MASS)
y <- rdiri(214, runif(3, 1, 3))
x <- as.matrix(fgl[, 2:9])
x <- x / rowSums(x)
mod <- tflr(y, x, x)
mod
}

