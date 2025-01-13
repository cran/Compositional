\name{Cross-validation for the alpha-SCLS model}
\alias{cv.ascls}
\title{
Cross-validation for the \eqn{alpha}-SCLS model
}

\description{
Cross-validation for the \eqn{alpha}-SCLS model.
}

\usage{
cv.ascls(y, x, a = seq(0.1, 1, by = 0.1), nfolds = 10, folds = NULL, seed = NULL)
}

\arguments{
\item{y}{
A numerical matrix with the simplicial response data. Zero values are allowed.
}
\item{x}{
A matrix with the simplicial predictor variables. Zero values are allowed.
}
\item{a}{
A vector or a single number of values of the \eqn{\alpha}-parameter. This has to be different from zero, and it can take negative values if there are no zeros in the simplicial response (y).
}
\item{nfolds}{
The number of folds for the K-fold cross validation, set to 10 by default.
}
\item{folds}{
If you have the list with the folds supply it here. You can also leave it NULL and it will create folds.
}
\item{seed}{
You can specify your own seed number here or leave it NULL.
}
}

\details{
The K-fold cross validation is performed in order to select the optimal value for \eqn{\alpha} of the \eqn{\alpha}-SCLS model.
}

\value{
A list including:
\item{runtime}{
The runtime of the cross-validation procedure.
}
\item{kl}{
The Kullback-Leibler divergence for every value of \eqn{\alpha}.
}
\item{js}{
The Jensen-Shannon divergence for every value of \eqn{\alpha}.
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
\code{\link{ascls}, \link{cv.atflr}
}
}

\examples{
library(MASS)
y <- rdiri( 214, runif(4, 1, 3) )
x <- as.matrix( fgl[, 2:9] )
mod <- cv.ascls(y, x, nfolds = 5)
}

