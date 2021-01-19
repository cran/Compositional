\name{Transformation-free linear regression for compositional responses and predictors}
\alias{tflr}
\title{
Transformation-free linear regression for compositional responses and predictors
}

\description{
Transformation-free linear regression for compositional responses and predictors.
}

\usage{
tflr(y, x, xnew = NULL)
}

\arguments{
\item{y}{
A matrix with the compositional response. Zero values are allowed.
}
\item{x}{
A matrix with the compositional predictors. Zero values are allowed.
}
\item{xnew}{
If you have new data use it, otherwise leave it NULL.
}
}

\details{
The transformation-free linear regression for compositional responses and predictors is implemented.
}

\value{
A list including:
\item{runtime}{
The time required by the regression.
}
\item{loglik}{
The log-likelihood \eqn{-\sum_{i=1}^ny_i\log{y_i/\hat{y}_i}}.
}
\item{be}{
The beta coefficients.
}
\item{est}{
The fitted values of xnew if xnew is not NULL.
}
}

\references{
Jacob Fiksel, Scott Zeger and Abhirup Datta (2020). A transformation-free linear regression for compositional outcomes and predictors.
https://arxiv.org/pdf/2004.07881.pdf
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{cv.tflr}, \link{ols.compcomp} \link{kl.alfapcr}
}
}

\examples{
library(MASS)
y <- rdiri(214, runif(3, 1, 3))
x <- as.matrix(fgl[, 2:9])
x <- x / rowSums(x)
mod <- tflr(y, x, x)
mod
}

