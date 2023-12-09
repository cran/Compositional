\name{Cross validation for the alpha-k-NN regression with compositional predictor variables}
\alias{alfaknnreg.tune}
\title{
Cross validation for the \eqn{\alpha}-k-NN regression with compositional predictor variables
}

\description{
Cross validation for the \eqn{\alpha}-k-NN regression with compositional predictor variables.
}

\usage{
alfaknnreg.tune(y, x, a = seq(-1, 1, by = 0.1), k = 2:10, nfolds = 10,
apostasi = "euclidean", method = "average", folds = NULL, seed = NULL, graph = FALSE)
}

\arguments{
\item{y}{
The response variable, a numerical vector.
}
\item{x}{
A matrix with the available compositional data. Zeros are allowed.
}
\item{a}{
A vector with a grid of values of the power transformation, it has to be between -1 and 1. If zero values are present it has to be greater than 0.
If \eqn{\alpha=0} the isometric log-ratio transformation is applied.
}
\item{k}{
The number of nearest neighbours to consider. It can be a single number or a vector.
}
\item{nfolds}{
The number of folds. Set to 10 by default.
}
\item{apostasi}{
The type of distance to use, either "euclidean" or "manhattan".
}
\item{method}{
If you want to take the average of the reponses of the k closest observations, type "average".
For the median, type "median" and for the harmonic mean, type "harmonic".
}
\item{folds}{
If you have the list with the folds supply it here. You can also leave it NULL and it will create folds.
}
\item{seed}{
If seed is TRUE the results will always be the same.
}
\item{graph}{
If graph is TRUE (default value) a filled contour plot will appear.
}
}

\details{
A k-fold cross validation for the \eqn{\alpha}-k-NN regression for compositional response data is performed.
}

\value{
A list including:
\item{mspe}{
The mean square error of prediction.
}
\item{performance}{
The minimum mean square error of prediction.
}
\item{opt_a}{
The optimal value of \eqn{\alpha}.
}
\item{opt_k}{
The optimal value of k.
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
\code{ \link{alfa.rda}, \link{alfa.fda} }
}

\examples{
library(MASS)
x <- as.matrix(fgl[, 2:9])
x <- x / rowSums(x)
y <- fgl[, 1]
mod <- alfaknnreg.tune(y, x, a = seq(0.2, 0.4, by = 0.1), k = 2:4, nfolds = 5)
}


