\name{Ridge regression with compositional data in the covariates side using the alpha-transformation}
\alias{alfa.ridge}
\title{
Ridge regression with compositional data in the covariates side using the \eqn{\alpha}-transformation
}

\description{
Ridge regression with compositional data in the covariates side using the \eqn{\alpha}-transformation.
}

\usage{
alfa.ridge(y, x, a, lambda, B = 1, xnew = NULL)
}

\arguments{
\item{y}{
A numerical vector containing the response variable values. If they are percentages, they are mapped onto \eqn{R} using the logit transformation.
}
\item{x}{
A matrix with the predictor variables, the compositional data. Zero values are allowed, but you must be carefull to choose strictly positive vcalues of \eqn{\alpha}.
}
\item{a}{
The value of the power transformation, it has to be between -1 and 1. If zero values are present it has to be greater than 0. If \eqn{\alpha=0} the isometric log-ratio transformation is applied.
}
\item{lambda}{
The value of the regularisation parameter, \eqn{\lambda}.
}
\item{B}{
If B > 1 bootstrap estimation of the standard errors is implemented.
}
\item{xnew}{
A matrix containing the new compositional data whose response is to be predicted. If you have no new data, leave this NULL as is by default.
}
}

\details{
The \eqn{\alpha}-transformation is applied to the compositional data first and then ridge components regression is performed.
}

\value{
The output of the \link{ridge.reg}.
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
\code{\link{ridge.reg}, \link{alfaridge.tune}, \link{alfaridge.plot}
}
}

\examples{
library(MASS)
y <- as.vector(fgl[, 1])
x <- as.matrix(fgl[, 2:9])
x <- x/ rowSums(x)
mod1 <- alfa.ridge(y, x, a = 0.5, lambda = 0.1, B = 1, xnew = NULL)
mod2 <- alfa.ridge(y, x, a = 0.5, lambda = 1, B = 1, xnew = NULL)
}

\keyword{ \eqn{\alpha}-transformation }
\keyword{ ridge regression }
