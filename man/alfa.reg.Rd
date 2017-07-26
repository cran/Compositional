\name{Regression with compositional data using the alpha-transformation}
\alias{alfa.reg}
\title{
Regression with compositional data using the \eqn{\alpha}-transformation
}

\description{
Regression with compositional data using the \eqn{\alpha}-transformation.
}

\usage{
alfa.reg(y, x, a, xnew = NULL, yb = NULL)
}

\arguments{
\item{y}{
A matrix with the compositional data.
}
\item{x}{
A matrix with the continuous predictor variables or a data frame including categorical predictor variables.
}
\item{a}{
The value of the power transformation, it has to be between -1 and 1. If zero values are present it has to be greater than 0. If \eqn{\alpha=0} the isometric log-ratio transformation is applied and the solution exists in a closed form, since it the classical mutivariate regression.
}
\item{xnew}{
If you have new data use it, otherwise leave it NULL.
}
\item{yb}{
If you have already transformed the data using the \eqn{\alpha}-transformation with the same \eqn{\alpha} as given in the argument "a", put it here. Othewrise leave it NULL.

This is intended to be used in the function \code{\link{alfareg.tune}} in order to speed up the process. The time difference in that function is small for small samples. But, if you have a few thousands and or a few more components, there will be bigger differences.
}
}

\details{
The \eqn{\alpha}-transformation is applied to the compositional data first and then multivariate regression is applied. This involves numerical optimisation.
}

\value{
A list including:
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
The fitted or the predicted values (if xnew is not NULL).
}
}

\references{
Tsagris M. (2015). Regression analysis with compositional data containing zero values. Chilean Journal of Statistics, 6(2): 47-57. http://arxiv.org/pdf/1508.01913v1.pdf

Tsagris M.T., Preston S. and Wood A.T.A. (2011). A data-based power transformation for compositional data.
In Proceedings of the 4th Compositional Data Analysis Workshop, Girona, Spain.
http://arxiv.org/pdf/1106.1451.pdf

Mardia K.V., Kent J.T., and Bibby J.M. (1979). Multivariate analysis. Academic press.

Aitchison J. (1986). The statistical analysis of compositional data. Chapman & Hall.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
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

\keyword{ \eqn{\alpha}-transformation }
\keyword{ multivariate regression }
