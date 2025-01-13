\name{Cross-validation for the SCRQ model}
\alias{cv.scrq}
\title{
Cross-validation for the SCRQ model
}

\description{
Cross-validation for the SCRQ model.
}

\usage{
cv.scrq(y, x, nfolds = 10, folds = NULL, seed = NULL)
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
You can specify your own seed number here or leave it NULL.
}
}

\details{
The function performs k-fold cross-validation for the absolute regression
where the beta coefficients are constained to be positive and sum to 1.
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
\code{\link{scrq}, \link{cv.scls}, \link{cv.tflr} }
}

\examples{
y <- rdiri(500, runif(3, 1, 3))
x <- rdiri(500, runif(3, 1, 3))
mod <- cv.scrq(y, x, nfolds = 5)
}
