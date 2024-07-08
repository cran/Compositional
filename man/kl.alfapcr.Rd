\name{Divergence based regression for compositional data with compositional data in the covariates side using the alpha-transformation}
\alias{kl.alfapcr}
\title{
Divergence based regression for compositional data with compositional data in the covariates side using the \eqn{\alpha}-transformation
}

\description{
Divergence based regression for compositional data with compositional data in the covariates side using the \eqn{\alpha}-transformation.
}

\usage{
kl.alfapcr(y, x, covar = NULL, a, k, xnew = NULL, B = 1, ncores = 1, tol = 1e-07,
maxiters = 50)
}

\arguments{
\item{y}{
A numerical matrixc with compositional data with or without zeros.
}
\item{x}{
A matrix with the predictor variables, the compositional data. Zero values are allowed.
}
\item{covar}{
If you have other covariates as well put themn here.
}
\item{a}{
The value of the power transformation, it has to be between -1 and 1. If zero values are present it has to be greater than 0.
If \eqn{\alpha=0} the isometric log-ratio transformation is applied.
}
\item{k}{
A number at least equal to 1. How many principal components to use.
}
\item{xnew}{
A matrix containing the new compositional data whose response is to be predicted. If you have no new data,
leave this NULL as is by default.
}
\item{B}{
If B is greater than 1 bootstrap estimates of the standard error are returned. If B=1, no standard errors are returned.
}
\item{ncores}{
If ncores is 2 or more parallel computing is performed. This is to be used for the case of bootstrap.
If B=1, this is not taken into consideration.
}
\item{tol}{
The tolerance value to terminate the Newton-Raphson procedure.
}
\item{maxiters}{
The maximum number of Newton-Raphson iterations.
}
}

\details{
The \eqn{\alpha}-transformation is applied to the compositional data first, the first k principal component scores are calcualted and used as predictor variables for the Kullback-Leibler divergence based regression model.
}

\value{
A list including:
\item{runtime}{
The time required by the regression.
}
\item{iters}{
The number of iterations required by the Newton-Raphson in the kl.compreg function.
}
\item{loglik}{
The log-likelihood. This is actually a quasi multinomial regression. This is bascially minus the half deviance, or
\eqn{- sum_{i=1}^ny_i\log{y_i/\hat{y}_i}}.
}
\item{be}{
The beta coefficients.
}
\item{seb}{
The standard error of the beta coefficients, if bootstrap is chosen, i.e. if B > 1.
}
\item{est}{
The fitted values of xnew if xnew is not NULL.
}
}

\references{
Alenazi A. (2019). Regression for compositional data with compositional data as predictor variables with or without zero values.
Journal of Data Science, 17(1): 219-238.
https://jds-online.org/journal/JDS/article/136/file/pdf

Tsagris M. (2015). Regression analysis with compositional data containing zero values. Chilean Journal of Statistics, 6(2): 47-57.
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
\code{\link{klalfapcr.tune}, \link{tflr}, \link{glm.pcr}, \link{alfapcr.tune}
}
}

\examples{
library(MASS)
y <- rdiri(214, runif(4, 1, 3))
x <- as.matrix(fgl[, 2:9])
x <- x / rowSums(x)
mod <- alfa.pcr(y = y, x = x, a = 0.7, k = 1)
mod
}

