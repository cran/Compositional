\name{Tuning of the divergence based regression for compositional data with compositional data in the covariates side using the alpha-transformation}
\alias{klalfapcr.tune}
\title{
Tuning of the divergence based regression for compositional data with compositional data in the covariates side using the \eqn{\alpha}-transformation
}

\description{
Tuning of the divergence based regression for compositional data with compositional data in the covariates side using the \eqn{\alpha}-transformation.
}

\usage{
klalfapcr.tune(y, x, covar = NULL, nfolds = 10, maxk = 50, a = seq(-1, 1, by = 0.1),
folds = NULL, graph = FALSE, tol = 1e-07, maxiters = 50, seed = NULL)
}

\arguments{
\item{y}{
A numerical matrix with compositional data with or without zeros.
}
\item{x}{
A matrix with the predictor variables, the compositional data. Zero values are allowed.
}
\item{covar}{
If you have other continuous covariates put themn here.
}
\item{nfolds}{
The number of folds for the K-fold cross validation, set to 10 by default.
}
\item{maxk}{
The maximum number of principal components to check.
}
\item{a}{
The value of the power transformation, it has to be between -1 and 1. If zero values are present it has to be greater than 0.
If \eqn{\alpha=0} the isometric log-ratio transformation is applied.
}
\item{folds}{
If you have the list with the folds supply it here. You can also leave it NULL and it will create folds.
}
\item{graph}{
If graph is TRUE (default value) a plot will appear.
}
\item{tol}{
The tolerance value to terminate the Newton-Raphson procedure.
}
\item{maxiters}{
The maximum number of Newton-Raphson iterations.
}
\item{seed}{
You can specify your own seed number here or leave it NULL.
}
}

\details{
The M-fold cross validation is performed in order to select the optimal values for \eqn{\alpha} and k, the number of principal components.
The \eqn{\alpha}-transformation is applied to the compositional data first, the first k principal component scores are calcualted and used as predictor variables for the Kullback-Leibler divergence based regression model. This procedure is performed M times during the M-fold cross validation.
}

\value{
A list including:
\item{mspe}{
A list with the KL divergence for each value of \eqn{\alpha} and k in every fold.
}
\item{performance}{
A matrix with the KL divergence for each value of \eqn{\alpha} averaged over all folds. If graph is set to TRUE this matrix is plotted.
}
\item{best.perf}{
The minimum KL divergence.
}
\item{params}{
The values of \eqn{\alpha} and k corresponding to the minimum KL divergence.
}
}

\references{
Alenazi A. (2019). Regression for compositional data with compositional data as predictor variables with or without zero values.
Journal of Data Science, 17(1): 219--238.
https://jds-online.org/journal/JDS/article/136/file/pdf

Tsagris M. (2015). Regression analysis with compositional data containing zero values. Chilean Journal of Statistics, 6(2): 47--57.
http://arxiv.org/pdf/1508.01913v1.pdf

Tsagris M.T., Preston S. and Wood A.T.A. (2011). A data-based power transformation for compositional data.
In Proceedings of the 4th Compositional Data Analysis Workshop, Girona, Spain.
http://arxiv.org/pdf/1106.1451.pdf
}

\author{
Initial code by Abdulaziz Alenazi. Modifications by Michail Tsagris.

R implementation and documentation: Abdulaziz Alenazi \email{a.alenazi@nbu.edu.sa}
and Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{kl.alfapcr}, \link{cv.tflr}, \link{glm.pcr}, \link{alfapcr.tune}
}
}

\examples{
library(MASS)
y <- rdiri( 214, runif(4, 1, 3) )
x <- as.matrix( fgl[, 2:9] )
x <- x / rowSums(x)
mod <- klalfapcr.tune(y = y, x = x, a = c(0.7, 0.8) )
mod
}

