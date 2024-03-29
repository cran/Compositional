\name{Multivariate linear regression}
\alias{multivreg}
\title{
Multivariate linear regression
}

\description{
Multivariate linear regression.
}

\usage{
multivreg(y, x, plot = TRUE, xnew = NULL)
}

\arguments{
\item{y}{
A matrix with the Eucldidean (continuous) data.
}
\item{x}{
A matrix with the predictor variable(s), they have to be continuous.
}
\item{plot}{
Should a plot appear or not?
}
\item{xnew}{
If you have new data use it, otherwise leave it NULL.
}
}

\details{
The classical multivariate linear regression model is obtained.
}

\value{
A list including:
\item{suma}{
A summary as produced by \code{\link{lm}}, which includes the coefficients, their standard error, t-values, p-values.
}
\item{r.squared}{
The value of the \eqn{R^2} for each univariate regression.
}
\item{resid.out}{
A vector with number indicating which vectors are potential residual outliers.
}
\item{x.leverage}{
A vector with number indicating which vectors are potential outliers in the predictor variables space.
}
\item{out}{
A vector with number indicating which vectors are potential outliers in the residuals and in the predictor variables space.
}
\item{est}{
The predicted values if xnew is not NULL.
}
}

\references{
K.V. Mardia, J.T. Kent and J.M. Bibby (1979). Multivariate Analysis. Academic Press.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} and Giorgos Athineou <gioathineou@gmail.com>.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{diri.reg}, \link{js.compreg}, \link{kl.compreg}, \link{ols.compreg}, \link{comp.reg}
}
}

\examples{
library(MASS)
x <- as.matrix(iris[, 1:2])
y <- as.matrix(iris[, 3:4])
multivreg(y, x, plot = TRUE)
}

\keyword{ multivariate linear regression }
