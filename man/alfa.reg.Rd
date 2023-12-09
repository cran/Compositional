\name{Regression with compositional data using the alpha-transformation}
\alias{alfa.reg}
\alias{alfa.reg2}
\alias{alfa.reg3}
\title{
Regression with compositional data using the \eqn{\alpha}-transformation
}

\description{
Regression with compositional data using the \eqn{\alpha}-transformation.
}

\usage{
alfa.reg(y, x, a, xnew = NULL, yb = NULL, seb = FALSE)
alfa.reg2(y, x, a, xnew = NULL, yb = NULL, seb = FALSE)
alfa.reg3(y, x, a = c(-1, 1), xnew = NULL, seb = FALSE)
}

\arguments{
\item{y}{
A matrix with the compositional data.
}
\item{x}{
A matrix with the continuous predictor variables or a data frame including categorical predictor variables.
}
\item{a}{
The value of the power transformation, it has to be between -1 and 1. If zero values are present it has to be greater than 0. If \eqn{\alpha=0} the isometric log-ratio transformation is applied and the solution exists in a closed form, since it the classical mutivariate regression. For the alfa.reg2() this should be a vector of \eqn{\alpha} values and the function call repeatedly the alfa.reg() function. For the alfa.reg3() function it should be a vector with two values, the endpoints of the interval of \eqn{\alpha}. This function searches for the optimal vaue of \eqn{\alpha} that minimizes the sum of squares of the errors. Using the \code{\link{optimize}} function it searches for the optimal value of \eqn{\alpha}. Instead of choosing the value of \eqn{\alpha} using \code{\link{alfareg.tune}} (that uses cross-validation) one can select it this way.
}
\item{xnew}{
If you have new data use it, otherwise leave it NULL.
}
\item{yb}{
If you have already transformed the data using the \eqn{\alpha}-transformation with the same \eqn{\alpha} as given in the argument "a", put it here. Othewrise leave it NULL.

This is intended to be used in the function \code{\link{alfareg.tune}} in order to speed up the process. The time difference in that function is small for small samples.
But, if you have a few thousands and or a few more components, there will be bigger differences.
}
\item{seb}{
Do you want the standard error of the coefficients to be returned? In the \code{\link{alfareg.tune}} function this extra computation that is avoided can save some time.
}
}

\details{
The \eqn{\alpha}-transformation is applied to the compositional data first and then multivariate regression is applied. This involves numerical optimisation.
}

\value{
For the alfa.reg() function a list including:
\item{runtime}{
The time required by the regression.
}
\item{be}{
The beta coefficients.
}
\item{seb}{
The standard error of the beta coefficients.
}
\item{est}{
The fitted values for xnew if xnew is not NULL.
}

For the alfa.reg2() function a list with as many sublists as the number of values of \eqn{\alpha}. Each element (sublist) of the list contains the above outcomes of the alfa.reg() function.

For the alfa.reg3() function a list with all previous elements plus an output "alfa", the optimal value of \eqn{\alpha}.
}

\references{
Tsagris M. (2015). Regression analysis with compositional data containing zero values. Chilean Journal of Statistics, 6(2): 47-57.
https://arxiv.org/pdf/1508.01913v1.pdf

Tsagris M.T., Preston S. and Wood A.T.A. (2011). A data-based power transformation for compositional data.
In Proceedings of the 4th Compositional Data Analysis Workshop, Girona, Spain.
https://arxiv.org/pdf/1106.1451.pdf

Mardia K.V., Kent J.T., and Bibby J.M. (1979). Multivariate analysis. Academic press.

Aitchison J. (1986). The statistical analysis of compositional data. Chapman & Hall.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} and Giorgos Athineou <gioathineou@gmail.com>.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{alfareg.tune}, \link{diri.reg}, \link{js.compreg}, \link{kl.compreg}, \link{ols.compreg}, \link{comp.reg}
}
}

\examples{
library(MASS)
x <- as.vector(fgl[1:40, 1])
y <- as.matrix(fgl[1:40, 2:9])
y <- y / rowSums(y)
mod <- alfa.reg(y, x, 0.2)
}

\keyword{ multivariate regression }
