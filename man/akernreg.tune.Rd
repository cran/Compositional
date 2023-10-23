\name{Cross validation for the alpha-kernel regression with compositional response data}
\alias{akernreg.tune}
\title{
Cross validation for the \eqn{\alpha}-kernel regression with compositional response data
}

\description{
Cross validation for the \eqn{\alpha}-kernel regression with compositional response data.
}

\usage{
akernreg.tune(y, x, a = seq(0.1, 1, by = 0.1), h = seq(0.1, 1, length = 10),
type = "gauss", nfolds = 10, folds = NULL, seed = NULL)
}

\arguments{
\item{y}{
A matrix with the compositional response data. Zeros are allowed.
}
\item{x}{
A matrix with the available predictor variables.
}
\item{a}{
A vector with a grid of values of the power transformation, it has to be between -1 and 1. If zero values are present it has to be greater than 0. If \eqn{\alpha=0} the isometric log-ratio transformation is applied.
}
\item{h}{
A vector with the bandwidth value(s) to consider.
}
\item{type}{
The type of kernel to use, "gauss" or "laplace".
}
\item{nfolds}{
The number of folds. Set to 10 by default.
}
\item{folds}{
If you have the list with the folds supply it here. You can also leave it NULL and it will create folds.
}
\item{seed}{
You can specify your own seed number here or leave it NULL.
}
}

\details{
A k-fold cross validation for the \eqn{\alpha}-kernel regression for compositional response data is performed.
}

\value{
A list including:
\item{kl}{
The Kullback-Leibler divergence for all combinations of \eqn{\alpha} and \eqn{h}.
}
\item{js}{
The Jensen-Shannon divergence for all combinations of \eqn{\alpha} and \eqn{h}.
}
\item{klmin}{
The minimum Kullback-Leibler divergence.
}
\item{jsmin}{
The minimum Jensen-Shannon divergence.
}
\item{kl.alpha}{
The optimal \eqn{\alpha} that leads to the minimum Kullback-Leibler divergence.
}
\item{kl.h}{
The optimal \eqn{h} that leads to the minimum Kullback-Leibler divergence.
}
\item{js.alpha}{
The optimal \eqn{\alpha} that leads to the minimum Jensen-Shannon divergence.
}
\item{js.h}{
The optimal \eqn{h} that leads to the minimum Jensen-Shannon divergence.
}
\item{runtime}{
The runtime of the cross-validation procedure.
}
}

\references{
Tsagris M., Alenazi A. and Stewart C. (2023).
Flexible non-parametric regression models for compositional response data with zeros.
Statistics and Computing, 33(5): 1--17.

https://link.springer.com/article/10.1007/s11222-023-10277-5
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{akern.reg}, \link{aknnreg.tune}, \link{aknn.reg}, \link{alfa.rda}, \link{alfa.fda}, \link{rda.tune} }
}

\examples{
y <- as.matrix( iris[, 1:3] )
y <- y / rowSums(y)
x <- iris[, 4]
mod <- akernreg.tune(y, x, a = c(0.4, 0.6), h = c(0.1, 0.2), nfolds = 5)
}

