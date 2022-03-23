\name{Cross validation for the transformation-free linear regression for compositional responses and predictors}
\alias{cv.tflr}
\title{
Cross validation for the transformation-free linear regression for compositional
responses and predictors
}

\description{
Cross validation for the transformation-free linear regression for compositional
responses and predictors.
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

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{tflr}, \link{cv.olscompcomp}, \link{klalfapcr.tune}
}
}

\examples{
\donttest{
library(MASS)
y <- rdiri(214, runif(3, 1, 3))
x <- as.matrix(fgl[, 2:9])
x <- x / rowSums(x)
mod <- cv.tflr(y, x)
mod
}
}

