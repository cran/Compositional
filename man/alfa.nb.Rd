\name{Naive Bayes classifiers for compositional data using the alpha-transformation}
\alias{alfa.nb}
\title{
Naive Bayes classifiers for compositional data using the \eqn{\alpha}-transformation
}

\description{
Naive Bayes classifiers for compositional data using the \eqn{\alpha}-transformation.
}

\usage{
alfa.nb(xnew, x, ina, a, type = "gaussian")
}
\arguments{
\item{xnew}{
A matrix with the new compositional predictor data whose class you want to predict. Zeros are allowed.
}
\item{x}{
A matrix with the available compositional predictor data. Zeros are allowed.
}
\item{ina}{
A vector of data. The response variable, which is categorical (factor is acceptable).
}
\item{a}{
This can be a vector of values or a single number.
}
\item{type}{
The type of naive Bayes, "gaussian", "cauchy" or "laplace".
}
}

\details{
The \eqn{\alpha}-transformation is applied to the compositional and a naive Bayes classifier is employed.
}

\value{
A matrix with the estimated groups. One column for each value of \eqn{\alpha}.
}

\references{
Tsagris M.T., Preston S. and Wood A.T.A. (2011). A data-based power transformation for compositional data.
In Proceedings of the 4th Compositional Data Analysis Workshop, Girona, Spain.
https://arxiv.org/pdf/1106.1451.pdf

Friedman J., Hastie T. and Tibshirani R. (2017). The elements of statistical learning. New York: Springer.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{comp.nb}, \link{alfa.rda}, \link{alfa.knn}, \link{comp.knn}, \link{mix.compnorm}
}
}

\examples{
x <- Compositional::rdiri(100, runif(5) )
ina <- rbinom(100, 1, 0.5) + 1
mod <- alfa.nb(x, x, a = c(0, 0.1, 0.2), ina )
}

