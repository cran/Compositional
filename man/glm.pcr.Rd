\name{Principal component generalised linear models}
\alias{pcr}
\alias{glm.pcr}
\title{
Principal component generalised linear models
}

\description{
Principal component generalised linear models.
}

\usage{
pcr(y, x, k = 1, xnew = NULL)
glm.pcr(y, x, k = 1, xnew = NULL)
}

\arguments{
\item{y}{
A numerical vector, a real values vector or a numeric vector with 0 and 1 (binary) or a vector with discrete (count) data.
}
\item{x}{
A matrixc with the predictor variable(s), they have to be continuous.
}
\item{k}{
A number greater than or equal to 1. How many principal components to use. In the case of "pcr" this can be a single number or a vector. In the second case you get results for the sequence of principal components.
}
\item{xnew}{
If you have new data use it, otherwise leave it NULL.
}
}

\details{
Principal component regression is performed with linear, binary logistic or Poisson regression, depending on the nature of the response variable. The principal components of the cross product of the independent variables are obtained and classical regression is performed. This is used in the function \code{\link{alfa.pcr}}.
}

\value{
A list including:
\item{be}{
The beta coefficients of the predictor variables computed via the principcal components if "pcr" is used.
}
\item{model}{
The summary of the logistic or Poisson regression model.
}
\item{per}{
The percentage of variance of the predictor variables retained by the k principal components.
}
\item{vec}{
The principal components, the loadings.
}
\item{est}{
The fitted or the predicted values (if xnew is not NULL).
}
}

\references{
Aguilera A.M., Escabias M. and Valderrama M.J. (2006). Using principal components for estimating logistic regression with high-dimensional multicollinear data.
Computational Statistics \& Data Analysis 50(8): 1905-1924.

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
\code{\link{alfa.pcr}, \link{alfapcr.tune}
}
}

\examples{
library(MASS)
x <- as.matrix(fgl[, 2:9])
y <- as.vector(fgl[, 1])
mod1 <- pcr(y, x, 1)
mod2 <- pcr(y, x, 2)
mod <- pcr(y, x, k = 1:4)  ## many results at once

x <- as.matrix(iris[, 1:4])
y<- rbinom(150, 1, 0.6)
mod<- glm.pcr(y, x, k = 1)
}

