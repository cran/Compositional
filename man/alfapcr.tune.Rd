\name{Tuning the number of PCs in the PCR with compositional data using the alpha-transformation}
\alias{alfapcr.tune}
\title{
Tuning the number of PCs in the PCR with compositional data using the \eqn{\alpha}-transformation
}

\description{
This is a cross-validation procedure to decide on the number of principal components when using regression with compositional data (as predictor variables) using the \eqn{\alpha}-transformation.
}

\usage{
alfapcr.tune(y, x, model = "gaussian", nfolds = 10, maxk = 50, a = seq(-1, 1, by = 0.1),
folds = NULL, ncores = 1, graph = TRUE, col.nu = 15, seed = FALSE)
}

\arguments{
\item{y}{
A vector with either continuous, binary or count data.
}
\item{x}{
A matrix with the predictor variables, the compositional data. Zero values are allowed.
}
\item{model}{
The type of regression model to fit. The possible values are "gaussian", "binomial" and "poisson".
}
\item{nfolds}{
The number of folds for the K-fold cross validation, set to 10 by default.
}
\item{maxk}{
The maximum number of principal components to check.
}
\item{a}{
The value of the power transformation, it has to be between -1 and 1. If zero values are present it has to be greater than 0. If \eqn{\alpha=0} the isometric log-ratio transformation is applied and the solution exists in a closed form, since it the classical mutivariate regression. The estimated bias correction via the (Tibshirani and Tibshirani (2009) criterion is applied.
}
\item{folds}{
If you have the list with the folds supply it here. You can also leave it NULL and it will create folds.
}
\item{ncores}{
How many cores to use. If you have heavy computations or do not want to wait for long time more than 1 core (if available) is suggested. It is advisable to use it if you have many observations and or many variables, otherwise it will slow down th process.
}
\item{graph}{
If graph is TRUE (default value) a filled contour plot will appear.
}
\item{col.nu}{
A number parameter for the filled contour plot, taken into account only if graph is TRUE.
}
\item{seed}{
If seed is TRUE the results will always be the same.
}
}

\details{
The \eqn{\alpha}-transformation is applied to the compositional data first and the function "pcr.tune" or "glmpcr.tune" is called.
}

\value{
If graph is TRUE a filled contour will appear.
A list including:
\item{mspe}{
The MSPE where rows correspond to the \eqn{\alpha} values and the columns to the number of principal components.
}
\item{best.par}{
The best pair of \eqn{\alpha} and number of principal components.
}
\item{performance}{
The minimum mean squared error of prediction.
}
\item{runtime}{
The time required by the cross-validation procedure.
}
}

\references{
Tsagris M. (2015). Regression analysis with compositional data containing zero values. Chilean Journal of
Statistics, 6(2): 47-57.
https://arxiv.org/pdf/1508.01913v1.pdf

Tsagris M.T., Preston S. and Wood A.T.A. (2011). A data-based power transformation for compositional data.
In Proceedings of the 4th Compositional Data Analysis Workshop, Girona, Spain.
https://arxiv.org/pdf/1106.1451.pdf

Jolliffe I.T. (2002). Principal Component Analysis.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} and Giorgos Athineou <gioathineou@gmail.com>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{alfa}, \link{profile}, \link{alfa.pcr}, \link{pcr.tune}, \link{glmpcr.tune}, \link{glm}
}
}

\examples{
library(MASS)
y <- as.vector(fgl[, 1])
x <- as.matrix(fgl[, 2:9])
x <- x/ rowSums(x)
mod <- alfapcr.tune(y, x, nfolds = 10, maxk = 50, a = seq(-1, 1, by = 0.1) )
}

\keyword{ principal components regression }
\keyword{ tuning }
