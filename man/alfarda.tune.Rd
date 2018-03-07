\name{Cross validation for the regularised discriminant analysis with compositional data using the alpha-transformation}
\alias{alfarda.tune}
\title{
Cross validation for the regularised discriminant analysis with compositional data using the \eqn{\alpha}-transformation
}

\description{
Cross validation for the regularised discriminant analysis with compositional data using the \eqn{\alpha}-transformation.
}

\usage{
alfarda.tune(x, ina, a = seq(-1, 1, by = 0.1), M = 10,
gam = seq(0, 1, by = 0.1), del = seq(0, 1, by = 0.1),
ncores = 1, mat = NULL)
}

\arguments{
\item{x}{
A matrix with the available compositional data. Zeros are allowed.
}
\item{ina}{
A group indicator variable for the avaiable data.
}
\item{a}{
A vector with a grid of values of the power transformation, it has to be between -1 and 1. If zero values are present it has to be greater than 0. 
If \eqn{\alpha=0} the isometric log-ratio transformation is applied.
}
\item{M}{
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
\item{mat}{
You can specify your own folds by giving a mat, where each column is a fold. Each column contains indices of the observations. You can also leave it 
NULL and it will create folds.
}
}

\details{
A k-fold cross validation is performed.
}

\value{
A list including:
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
}

\references{
Friedman Jerome, Trevor Hastie and Robert Tibshirani (2009). The elements of statistical learning, 2nd edition. Springer, Berlin

Tsagris M.T., Preston S. and Wood A.T.A. (2016). Improved classification for compositional data using the \eqn{\alpha}-transformation. 

}

\author{
Michail Tsagris

R implementation and documentation: Giorgos Athineou <gioathineou@gmail.com> and Michail Tsagris <mtsagris@yahoo.gr>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{alfa.rda}, \link{rda.tune} }
}

\examples{
library(MASS)
x <- as.matrix(fgl[, 2:9])
x <- x / rowSums(x)
ina <- fgl[, 10]
moda <- alfarda.tune(x, ina, a = seq(0.7, 1, by = 0.1), M = 10,
gam = seq(0.1, 0.3, by = 0.1), del = seq(0.1, 0.3, by = 0.1),
ncores = 1, mat = NULL)
}

\keyword{ Regularised discriminant analysis  }
\keyword{ cross-validation }
