\name{Permutation linear independence test in the TFLR model}
\alias{tflr.indeptest}
\title{
Permutation linear independence test in the TFLR model
}

\description{
Permutation linear independence test in the TFLR model.
}

\usage{
tflr.indeptest(y, x, R = 999, ncores = 1)
}

\arguments{
\item{y}{
A matrix with the compositional data (dependent variable). Zero values are allowed.
}
\item{x}{
A matrix with the compositional predictors. Zero values are in general allowed, but there can be cases when these are problematic.
}
\item{R}{
The number of permutations to perform.
}
\item{ncores}{
The number of cores to use in case you are interested for parallel computations.
}
}

\details{
Permutation independence test in the constrained linear least squares for compositional
responses and predictors is performed. The observed test statistic is the Kullback-Leibler divergence computed by \code{\link{tflr}}. Then, the rows of X are permuted B times and each time the TFLR is performed and the Kullback-Leibler is computed. The p-value is then computed in the usual way.
}

\value{
The p-value for the test of linear independence between the simplicial response Y and the simplicial predictor X.
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
\code{\link{tflr}, \link{scls}, \link{tflr.betest} }
}

\examples{
library(MASS)
set.seed(1234)
y <- rdiri(214, runif(4, 1, 3))
x <- as.matrix(fgl[, 2:9])
x <- x / rowSums(x)
tflr.indeptest(y, x, R = 10)
}

