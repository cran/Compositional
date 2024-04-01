\name{Cross validation for the TFLR model}
\alias{cv.tflr}
\title{
Cross validation for the TFLR model
}

\description{
Cross validation for the TFLR model.
}

\usage{
cv.tflr(y, x, nfolds = 10, folds = NULL, seed = NULL)
}

\arguments{
\item{y}{
A matrix with compositional response data. Zero values are allowed.
}
\item{x}{
A matrix with compositional predictors. Zero values are allowed.
}
\item{nfolds}{
The number of folds to be used. This is taken into consideration only if the
folds argument is not supplied.
}
\item{folds}{
If you have the list with the folds supply it here. You can also leave it NULL
and it will create folds.
}
\item{seed}{
If seed is TRUE the results will always be the same.
}
}

\details{
A k-fold cross validation for the transformation-free linear regression for
compositional responses and predictors is performed.
}

\value{
A list including:
\item{runtime}{
The runtime of the cross-validation procedure.
}
\item{kl}{
The Kullback-Leibler divergences for all runs.
}
\item{js}{
The Jensen-Shannon divergences for all runs.
}
\item{perf}{
The average Kullback-Leibler divergence and average Jensen-Shannon divergence.
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
\code{ \link{tflr}, \link{cv.scls}, \link{klalfapcr.tune}
}
}

\examples{
\donttest{
library(MASS)
y <- rdiri(100, runif(3, 1, 3))
x <- as.matrix(fgl[1:100, 2:9])
x <- x / rowSums(x)
mod <- cv.tflr(y, x)
mod
}
}

