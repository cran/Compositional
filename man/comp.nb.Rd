\name{Naive Bayes classifiers for compositional data}
\alias{comp.nb}
\title{
Naive Bayes classifiers for compositional data
}

\description{
Naive Bayes classifiers for compositional data.
}

\usage{
comp.nb(xnew = NULL, x, ina, type = "beta")
}
\arguments{
\item{xnew}{
A matrix with the new compositional predictor data whose class you want to predict. Zeros are not allowed
}
\item{x}{
A matrix with the available compositional predictor data. Zeros are not allowed
}
\item{ina}{
A vector of data. The response variable, which is categorical (factor is acceptable).
}
\item{type}{
The type of naive Bayes, "beta", "logitnorm", "cauchy", "laplace", "gamma", "normlog" or "weibull". For the last 4 distributions, the negative of the logarithm of the compositional data is applied first.
}
}

%\details{
%
%}

\value{
Depending on the classifier a list including (the ni and est are common for all classifiers):
\item{shape}{
A matrix with the shape parameters.
}
\item{scale}{
A matrix with the scale parameters.
}
\item{expmu}{
A matrix with the mean parameters.
}
\item{sigma}{
A matrix with the (MLE, hence biased) variance parameters.
}
\item{location}{
A matrix with the location parameters (medians).
}
\item{scale}{
A matrix with the scale parameters.
}
\item{mean}{
A matrix with the scale parameters.
}
\item{var}{
A matrix with the variance parameters.
}
\item{a}{
A matrix with the "alpha" parameters.
}
\item{b}{
A matrix with the "beta" parameters.
}
\item{ni}{
The sample size of each group in the dataset.
}
\item{est}{
The estimated group of the xnew observations. It returns a numerical value back regardless of the target
variable being numerical as well or factor. Hence, it is suggested that you do \"as.numeric(ina)\" in order to
see what is the predicted class of the new data.
}
}

\references{
Friedman J., Hastie T. and Tibshirani R. (2017). The elements of statistical learning.
New York: Springer.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{cv.compnb}, \link{alfa.rda}, \link{alfa.knn}, \link{comp.knn}, \link{mix.compnorm}, \link{dda}
}
}

\examples{
x <- Compositional::rdiri(100, runif(5) )
ina <- rbinom(100, 1, 0.5) + 1
a <- comp.nb(x, x, ina, type = "beta")
}

