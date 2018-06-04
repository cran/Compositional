\name{Multivariate or univariate regression with compositional data in the covariates side using the alpha-transformation}
\alias{alfa.pcr}
\title{
Multivariate or univariate regression with compositional data in the covariates side using the \eqn{\alpha}-transformation
}

\description{
Multivariate or univariate regression with compositional data in the covariates side using the \eqn{\alpha}-transformation.
}

\usage{
alfa.pcr(y, x, a, k, model = "gaussian", xnew = NULL)
}

\arguments{
\item{y}{
A numerical vector containing the response variable values. They can be continuous, binary, discrete (counts).
This can also be a vector with discrete values or a factor for the multinomial regression (model = "multinomial").
}
\item{x}{
A matrix with the predictor variables, the compositional data.
}
\item{a}{
The value of the power transformation, it has to be between -1 and 1. If zero values are present it has to be greater than 0.
If \eqn{\alpha=0} the isometric log-ratio transformation is applied.
}
\item{k}{
A number at least equal to 1. How many principal components to use.
}
\item{model}{
The type of regression model to fit. The possible values are "gaussian", "multinomial", "binomial" and "poisson".
}
\item{xnew}{
A matrix containing the new compositional data whose response is to be predicted. If you have no new data, leave this NULL as is by default.
}
}

\details{
The \eqn{\alpha}-transformation is applied to the compositional data first ,the first k principal component scores are calcualted and used as predictor variables for a regression model. The family of distributions can be either, "normal" for continuous response and hence normal distribution, "binomial" corresponding to binary response and hence logistic regression or "poisson" for count response and poisson regression.
}

\value{
A lis tincluding:
\item{be}{
If linear regression was fitted, the regression coefficients of the k principal component scores on the response variable y.
}
\item{mod}{
If another regression model was fitted its outcome as produced in the package \bold{Rfast}.
}
\item{per}{
The percentage of variance explained by the first k principal components.
}
\item{vec}{
The first k principal components, loadings or eigenvectors. These are useful for future prediction in the sense that one needs not fit the whole model again.
}
\item{est}{
If the argument "xnew" was given these are the predicted or estimated values, otherwise it is NULL.
}
}

\references{
Tsagris M. (2015). Regression analysis with compositional data containing zero values. Chilean Journal of Statistics, 6(2): 47-57.
http://arxiv.org/pdf/1508.01913v1.pdf

Tsagris M.T., Preston S. and Wood A.T.A. (2011). A data-based power transformation for compositional data.
In Proceedings of the 4th Compositional Data Analysis Workshop, Girona, Spain.
http://arxiv.org/pdf/1106.1451.pdf
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <gioathineou@gmail.com>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{pcr}, \link{glm.pcr}, \link{alfapcr.tune}
}
}

\examples{
library(MASS)
y <- as.vector(fgl[, 1])
x <- as.matrix(fgl[, 2:9])
x <- x / rowSums(x)
mod <- alfa.pcr(y = y, x = x, 0.7, 1)
mod
}

\keyword{ \eqn{\alpha}-transformation }
\keyword{ principal components regression }
