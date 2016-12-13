\name{Spatial median regression}
\alias{spatmed.reg}
\title{
Spatial median regression
}

\description{
Spatial median regression with Euclidean data.
}

\usage{
spatmed.reg(y, x, xnew = NULL, tol = 1e-07, ses = FALSE)
}

\arguments{
\item{y}{
A matrix with the compositional data. Zero values are not allowed.
}
\item{x}{
The predictor variable(s), they have to be continuous.
}
\item{xnew}{
If you have new data use it, otherwise leave it NULL.
}
\item{tol}{
The threshold upon which to stop the iterations of the Newton-Rapshon algorithm.
}
\item{ses}{
If you want to extract the standard errors of the parameters, set this to TRUE. Be careful though as this can slow down the algorithm dramatically.
In a run example with 10,000 observations and 10 variables for y and 30 for x, whe ses = FALSE the algorithm can take 0.20 seconds, but when ses = TRUE it can go up to 140 seconds.
}
}

\details{
The objective function is the minimization of the sum of the absolute residuals. It is the multivariate generalisation of the median regression.
This function is used by \code{\link{comp.reg}}.
}

\value{
A list including:
\item{iter}{
The number of iterations that were required.
}
\item{runtime}{
The time required by the regression.
}
\item{be}{
The beta coefficients.
}
\item{seb}{
The standard error of the beta coefficients is returned if ses=TRUE and NULL otherwise.
}
\item{est}{
The fitted or the predicted values (if xnew is not NULL).
}
}

\references{
Biman Chakraborty (2003) On multivariate quantile regression. Journal of Statistical Planning and Inference http://www.stat.nus.edu.sg/export/sites/dsap/research/documents/tr01_2000.pdf
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{multivreg}, \link{comp.reg}, \link{alfa.reg}, \link{js.compreg}, \link{diri.reg}
}
}

\examples{
library(MASS)
x <- as.matrix(iris[, 3:4])
y <- as.matrix(iris[, 1:2])
mod1 <- spatmed.reg(y, x)
mod2 <- multivreg(y, x, plot = FALSE)
}

\keyword{ spatial median regression }
\keyword{ multivariate regression }
