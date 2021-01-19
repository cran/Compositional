\name{Cross-validation for the constrained linear least squares for compositional responses and predictors}
\alias{cv.olscompcomp}
\title{
Cross-validation for the constrained linear least squares for compositional responses and predictors
}

\description{
Cross-validation for the constrained linear least squares for compositional responses and predictors.
}

\usage{
cv.olscompcomp(y, x, rs = 5, tol = 1e-4, nfolds = 10, folds = NULL, seed = FALSE)
}

\arguments{
\item{y}{
A matrix with compositional response data. Zero values are allowed.
}
\item{x}{
A matrix with compositional predictors. Zero values are allowed.
}
\item{rs}{
The number of times to run the constrained optimisation using different random starting values each time.
}
\item{tol}{
The threshold upon which to stop the iterations of the constrained optimisation.
}
\item{nfolds}{
The number of folds to be used. This is taken into consideration only if the folds argument is not supplied.
}
\item{folds}{
If you have the list with the folds supply it here. You can also leave it NULL and it will create folds.
}
\item{seed}{
If seed is TRUE the results will always be the same.
}
}

\details{
The function performs k-fold cross-validation for the least squares regression where the beta coefficients are constained to be positive and sum to 1.
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

%\references{
  %Murteira, Jose MR, and Joaquim JS Ramalho 2016. Regression analysis of multivariate fractional data. Econometric Reviews 35(4): 515-552.
  %}
\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
  %%  ~~further notes~~
    %}

\seealso{
\code{\link{ols.compcomp}, \link{cv.tflr}, \link{klalfapcr.tune} }
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

