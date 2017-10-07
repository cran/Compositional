\name{Principal components regression}
\alias{pcr}
\title{
Principal components regression
}

\description{
Principal components regression.
}

\usage{
pcr(y, x, k = 1, xnew = NULL)
}

\arguments{
  \item{y}{
  A real values vector.
}
  \item{x}{
  A matrix with the predictor variable(s), they have to be continuous.
}
  \item{k}{
The number of principal components to use.
}
  \item{xnew}{
  If you have new data use it, otherwise leave it NULL.
}
}

\details{
The principal components of the cross product of the independent variables are obtained and classical regression is performed. This is used in the function
\code{\link{alfa.pcr}}.
}

\value{
A list including:
  \item{beta}{
The beta coefficients.
}
  \item{parameters}{
The beta coefficients and their standard eror.
}
  \item{mse}{
The MSE of the linear regression, if xnew is NULL, becuase it needs the fitted values.
}
  \item{adj.rsq}{
The value of the adusted \eqn{R^2} if xnew is NULL.
}
  \item{per}{
The percentage of variance of the cross product of the independent variables explained by the k components.
}
  \item{est}{
The fitted or the predicted values (if xnew is not NULL).
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
\code{\link{pcr.tune}, \link{alfa.pcr}, \link{glm.pcr}
}
}

\examples{
library(MASS)
x <- as.matrix(fgl[, 2:9])
y <- as.vector(fgl[, 1])
mod1 <- pcr(y, x, 1)
mod2 <- pcr(y, x, 2)
}

\keyword{ principal components }
\keyword{ linear regression }
