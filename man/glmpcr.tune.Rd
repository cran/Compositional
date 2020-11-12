\name{Tuning the principal components with GLMs}
\alias{pcr.tune}
\alias{glmpcr.tune}
\alias{multinompcr.tune}
\title{
Tuning the principal components with GLMs
}

\description{
Tuning the number of principal components in the generalised linear models.
}

\usage{
pcr.tune(y, x, nfolds = 10, maxk = 50, folds = NULL, ncores = 1,
seed = FALSE, graph = TRUE)

glmpcr.tune(y, x, nfolds = 10, maxk = 10, folds = NULL, ncores = 1,
seed = FALSE, graph = TRUE)

multinompcr.tune(y, x, nfolds = 10, maxk = 10, folds = NULL, ncores = 1,
seed = FALSE, graph = TRUE)
}

\arguments{
\item{y}{
A real valued vector for "pcr.tune". A real valued vector for the "glmpcr.tune" with either two numbers, 0 and 1 for example, for the binomial regression or with positive discrete numbers for the poisson. For the "multinompcr.tune" a vector or a factor with more than just two values. This is a multinomial regression.
}
\item{x}{
A matrix with the predictor variables, they have to be continuous.
}
\item{nfolds}{
The number of folds in the cross validation.
}
\item{maxk}{
The maximum number of principal components to check.
}
\item{folds}{
If you have the list with the folds supply it here. You can also leave it NULL and it will create folds.
}
\item{ncores}{
The number of cores to use. If more than 1, parallel computing will take place. It is advisable to use it if you have many observations and or many variables, otherwise it will slow down th process.
}
\item{seed}{
If seed is TRUE the results will always be the same.
}
\item{graph}{
If graph is TRUE a plot of the performance for each fold along the values of \eqn{\alpha} will appear.
}
}

\details{
Cross validation is performed to select the optimal number of principal components in the GLMs or the multinomial regression. This is used
by \code{\link{alfapcr.tune}}.
}

\value{
If graph is TRUE a plot of the performance versus the number of principal components will appear.
A list including:
\item{msp}{
A matrix with the mean deviance of prediction or mean accuracy for every fold.
}
\item{mpd}{
A vector with the mean deviance of prediction or mean accuracy, each value corresponds to a number of principal components.
}
\item{k}{
The number of principal components which minimizes the deviance or maximises the accuracy.
}
\item{performance}{
The optimal performance, MSE for the linea regression, minimum deviance for the GLMs and maximum accuracy for the multinomial regression.
}
\item{runtime}{
The time required by the cross-validation procedure.
}
}

\references{
Aguilera A.M., Escabias M. and Valderrama M.J. (2006). Using principal components for estimating logistic regression with high-dimensional multicollinear data. Computational Statistics \& Data Analysis 50(8): 1905-1924.

Jolliffe I.T. (2002). Principal Component Analysis.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{pcr.tune}, \link{glm.pcr}, \link{alfa.pcr}, \link{alfapcr.tune}
}
}

\examples{
library(MASS)
x <- as.matrix(fgl[, 2:9])
y <- rpois(214, 10)
glmpcr.tune(y, x, nfolds = 10, maxk = 20, folds = NULL, ncores = 1)
}

\keyword{ GLM }
\keyword{ principal components }
