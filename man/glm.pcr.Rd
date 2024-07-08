\name{Principal component generalised linear models}
\alias{glm.pcr}
\title{
Principal component generalised linear models
}

\description{
Principal component generalised linear models.
}

\usage{
glm.pcr(y, x, k = 1, xnew = NULL)
}

\arguments{
\item{y}{
A numerical vector with 0 and 1 (binary) or a
vector with discrete (count) data.
}
\item{x}{
A matrix with the predictor variable(s), they have to be continuous.
}
\item{k}{
A number greater than or equal to 1. How many principal components to use. 
You may get results for the sequence of principal components.
}
\item{xnew}{
If you have new data use it, otherwise leave it NULL.
}
}

\details{
Principal component regression is performed with binary logistic or Poisson regression,
depending on the nature of the response variable. The principal components of the cross product
of the independent variables are obtained and classical regression is performed. This is used
in the function \code{\link{alfa.pcr}}.
}

\value{
A list including:
\item{model}{
The summary of the logistic or Poisson regression model as returned by the package Rfast.
}
\item{per}{
The percentage of variance of the predictor variables retained by the k principal components.
}
\item{vec}{
The principal components, the loadings.
}
\item{est}{
The fitted or the predicted values (if xnew is not NULL). If the argument \code{k} is a vector then this is a matrix with the estimated values for each number of components.
}
}

\references{
Aguilera A.M., Escabias M. and Valderrama M.J. (2006). Using principal components for estimating
logistic regression with high-dimensional multicollinear data.
Computational Statistics & Data Analysis 50(8): 1905-1924.

Jolliffe I.T. (2002). Principal Component Analysis.
}

\author{
Michail Tsagris.

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
x <- as.matrix(iris[, 1:4])
y <- rbinom(150, 1, 0.6)
mod <- glm.pcr(y, x, k = 1)
}

