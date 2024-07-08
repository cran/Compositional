\name{Cross validation for the alpha-k-NN regression with compositional response data}
\alias{aknnreg.tune}
\title{
Cross validation for the \eqn{\alpha}-k-NN regression with compositional response data
}

\description{
Cross validation for the \eqn{\alpha}-k-NN regression with compositional response data.
}

\usage{
aknnreg.tune(y, x, a = seq(0.1, 1, by = 0.1), k = 2:10, apostasi = "euclidean",
nfolds = 10, folds = NULL, seed = NULL, rann = FALSE)
}

\arguments{
\item{y}{
A matrix with the compositional response data. Zeros are allowed.
}
\item{x}{
A matrix with the available predictor variables.
}
\item{a}{
A vector with a grid of values of the power transformation, it has to be between -1 and 1. If zero values are present it has to be greater than 0.
If \eqn{\alpha=0} the isometric log-ratio transformation is applied.
}
\item{k}{
The number of nearest neighbours to consider. It can be a single number or a vector.
}
\item{apostasi}{
The type of distance to use, either "euclidean" or "manhattan".
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
\item{rann}{
If you have large scale datasets and want a faster k-NN search, you can use kd-trees implemented in the R package "Rnanoflann".
In this case you must set this argument equal to TRUE. Note however, that in this case, the only available distance is by default "euclidean".
}
}

\details{
A k-fold cross validation for the \eqn{\alpha}-k-NN regression for compositional response data is performed.
}

\value{
A list including:
\item{kl}{
The Kullback-Leibler divergence for all combinations of \eqn{\alpha} and \eqn{k}.
}
\item{js}{
The Jensen-Shannon divergence for all combinations of \eqn{\alpha} and \eqn{k}.
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
\item{kl.k}{
The optimal \eqn{k} that leads to the minimum Kullback-Leibler divergence.
}
\item{js.alpha}{
The optimal \eqn{\alpha} that leads to the minimum Jensen-Shannon divergence.
}
\item{js.k}{
The optimal \eqn{k} that leads to the minimum Jensen-Shannon divergence.
}
\item{runtime}{
The runtime of the cross-validation procedure.
}
}

\references{
Tsagris M., Alenazi A. and Stewart C. (2023).
Flexible non-parametric regression models for compositional response data with zeros.
Statistics and Computing, 33(106).

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
\code{\link{aknn.reg}, \link{akernreg.tune}, \link{akern.reg}, \link{alfa.rda}, \link{alfa.fda} }
}

\examples{
y <- as.matrix( iris[, 1:3] )
y <- y / rowSums(y)
x <- iris[, 4]
mod <- aknnreg.tune(y, x, a = c(0.4, 0.6), k = 2:4, nfolds = 5)
}

