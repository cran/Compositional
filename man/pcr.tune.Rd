\name{Tuning of the principal components regression}
\alias{pcr.tune}
\title{
Tuning of the principal components regression
}

\description{
Tuning the number of principal components in the principal components regression.
}

\usage{
pcr.tune(y, x, M = 10, maxk = 50, mat = NULL, ncores = 1, graph = TRUE)
}

\arguments{
\item{y}{
A real valued vector.
}
\item{x}{
A matrix with the predictor variables, they have to be continuous.
}
\item{M}{
The number of folds in the cross validation.
}
\item{maxk}{
The maximum number of principal components to check.
}
\item{mat}{
You can specify your own folds by giving a mat, where each column is a fold. Each column contains indices of the observations. 
You can also leave it NULL and it will create folds.
}
\item{ncores}{
The number of cores to use. If more than 1, parallel computing will take place. It is advisable to use it if you have many 
observations and or many variables, otherwise it will slow down th process.
}
\item{graph}{
If graph is TRUE a plot of the performance for each fold along the values of \eqn{\alpha} will appear.
}
}

\details{
Cross validation is performed to select the optimal number of principal components in the regression. 
This is used by \code{\link{alfapcr.tune}}.
}

\value{
A list including:
If graph is TRUE a plot of the performance versus the number of principal components will appear.
\item{msp}{
A matrix with the mean squared error of prediction (MSPE) for every fold.
}
\item{mspe}{
A vector with the mean squared error of prediction (MSPE), each value corresponds to a number of principal components.
}
\item{k}{
The number of principal components which minimizes the MSPE.
}
\item{performance}{
The lowest value of the MSPE.
}
\item{runtime}{
The time required by the cross-validation procedure.
}
}

\references{
Jolliffe I.T. (2002). Principal Component Analysis.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <gioathineou@gmail.com>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{glmpcr.tune}, \link{glm.pcr}, \link{alfa.pcr}, \link{alfapcr.tune}
}
}

\examples{
library(MASS)
x <- as.matrix(fgl[, 2:9])
y <- as.vector(fgl[, 1])
pcr.tune(y, x, M = 10, maxk = 50, mat = NULL, ncores = 1)
}

\keyword{ principal components }
\keyword{ linear regression }
