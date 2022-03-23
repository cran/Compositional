\name{Cross validation for the regularised and flexible discriminant analysis with compositional data using the alpha-transformation}
\alias{alfarda.tune}
\alias{alfafda.tune}
\title{
Cross validation for the regularised and flexible discriminant analysis with compositional data using the \eqn{\alpha}-transformation
}

\description{
Cross validation for the regularised and flexible discriminant analysis with compositional data using the \eqn{\alpha}-transformation.
}

\usage{
alfarda.tune(x, ina, a = seq(-1, 1, by = 0.1), nfolds = 10,
gam = seq(0, 1, by = 0.1), del = seq(0, 1, by = 0.1),
ncores = 1, folds = NULL, stratified = TRUE, seed = NULL)

alfafda.tune(x, ina, a = seq(-1, 1, by = 0.1), nfolds = 10,
folds = NULL, stratified = TRUE, seed = NULL, graph = FALSE)
}

\arguments{
\item{x}{
A matrix with the available compositional data. Zeros are allowed.
}
\item{ina}{
A group indicator variable for the avaiable data.
}
\item{a}{
A vector with a grid of values of the power transformation, it has to be between -1 and 1. If zero values are present it has to be greater than 0. If \eqn{\alpha=0} the isometric log-ratio transformation is applied.
}
\item{nfolds}{
The number of folds. Set to 10 by default.
}
\item{gam}{
A vector of values between 0 and 1. It is the weight of the pooled covariance and the diagonal matrix.
}
\item{del}{
A vector of values between 0 and 1. It is the weight of the LDA and QDA.
}
\item{ncores}{
The number of cores to use. If it is more than 1 parallel computing is performed. It is advisable to use it if you have many observations and or
many variables, otherwise it will slow down th process.
}
\item{folds}{
If you have the list with the folds supply it here. You can also leave it NULL and it will create folds.
}
\item{stratified}{
Do you want the folds to be created in a stratified way? TRUE or FALSE.
}
\item{seed}{
You can specify your own seed number here or leave it NULL.
}
\item{graph}{
If graph is TRUE (default value) a filled contour plot will appear.
}
}

\details{
A k-fold cross validation is performed.
}

\value{
For the alfa.rda a list including:
\item{res}{
The estimated optimal rate and the best values of \eqn{\alpha}, \eqn{gamma} and \eqn{delta}.
}
\item{percent}{
For the best value of \eqn{\alpha} the averaged over all folds best prates of correct classification. It is a matrix, where rows correspond to
the \eqn{\gamma} values and columns correspond to \eqn{\delta} values.
}
\item{se}{
The estimated standard errors of the "percent" matrix.
}
\item{runtime}{
The runtime of the cross-validation procedure.
}

For the alfa.fda a list including:
\item{per}{
The performance of the fda in each fold for each value of \eqn{\alpha}.
}
\item{performance}{
The average performance for each value of \eqn{\alpha}.
}
\item{opt_a}{
The optimal value of \eqn{\alpha}.
}
\item{runtime}{
The runtime of the cross-validation procedure.
}
}

\references{
Friedman Jerome, Trevor Hastie and Robert Tibshirani (2009). The elements of statistical learning, 2nd edition. Springer, Berlin

Tsagris M.T., Preston S. and Wood A.T.A. (2016). Improved classification for compositional data using the \eqn{\alpha}-transformation. Jounal of Classification, 33(2):243-261.

Hastie, Tibshirani and Buja (1994). Flexible Disriminant Analysis by Optimal Scoring. Journal of the American Statistical Association, 89(428):1255-1270.
}

\author{
Michail Tsagris.

R implementation and documentation: Giorgos Athineou <gioathineou@gmail.com> and 
Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{alfa.rda}, \link{alfanb.tune}, \link{cv.dda}, \link{compknn.tune}, \link{rda.tune}, \link{cv.compnb} }
}

\examples{
library(MASS)
x <- as.matrix(fgl[, 2:9])
x <- x / rowSums(x)
ina <- fgl[, 10]
moda <- alfarda.tune(x, ina, a = seq(0.7, 1, by = 0.1), nfolds = 10,
gam = seq(0.1, 0.3, by = 0.1), del = seq(0.1, 0.3, by = 0.1) )
}

\keyword{ Regularised discriminant analysis  }
\keyword{ cross-validation }
