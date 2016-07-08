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
  A numerical vector, either 0 and 1 (binary) or discrete (count) data.
}
  \item{x}{
  The predictor variable(s), they have to be continuous.
}
  \item{k}{
    A number at least equal to 1. How many principal components to use.
}
  \item{xnew}{
  If you have new data use it, otherwise leave it NULL.
}
}

\details{
Principal component regression is performed with binary logistic or Poisson regression, depending on the nature of the response variable.
}

\value{
A list including:
  \item{model}{
The summary of the glm model
}
  \item{per}{
The percentage of variance of the predictor variables retained by the k principal components.
}
  \item{est}{
The fitted or the predicted values (if xnew is not NULL).
}
}

\references{
Aguilera A.M., Escabias M. and Valderrama M.J. (2006). Using principal components for estimating logistic regression with high-dimensional multicollinear data.
Computational Statistics & Data Analysis 50(8): 1905-1924.

Jolliffe I.T. (2002). Principal Component Analysis.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{pcr}, \link{alfa.pcr}, \link{alfapcr.tune}, \link{glm}
}
}

\examples{
library(MASS)
x <- iris[, 1:4]
y<- rbinom(150, 1, 0.6)
mod<- glm.pcr(y, x, k = 1)
}

\keyword{ principal components regression }
\keyword{ generalised linear models }
