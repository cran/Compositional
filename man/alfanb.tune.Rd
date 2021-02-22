\name{Cross-validation for the naive Bayes classifiers for compositional data using the alpha-transformation}
\alias{alfanb.tune}
\title{
Cross-validation for the naive Bayes classifiers for compositional data using the \eqn{\alpha}-transformation
}

\description{
Cross-validation for the naive Bayes classifiers for compositional data using the \eqn{\alpha}-transformation.
}

\usage{
alfanb.tune(x, ina, a = seq(-1, 1, by = 0.1), type = "gaussian",
folds = NULL, nfolds = 10, stratified = TRUE, seed = FALSE)
}
\arguments{
\item{x}{
A matrix with the available data, the predictor variables.
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
\item{folds}{
A list with the indices of the folds.
}
\item{nfolds}{
The number of folds to be used. This is taken into consideration only if "folds" is NULL.
}
\item{stratified}{
Do you want the folds to be selected using stratified random sampling? This preserves the analogy of the samples of each group.
Make this TRUE if you wish.
}
\item{seed}{
If you set this to TRUE, the same folds will be created every time.
}
}

\details{
This function estimates the performance of the naive Bayes classifier for each value of \eqn{\alpha} of the \eqn{\alpha}-transformation.
}

\value{
A list including:
\item{crit}{
A vector whose length is equal to the number of k and is the accuracy metric for each k.
For the classification case it is the percentage of correct classification.
}
}

\references{
Friedman J., Hastie T. and Tibshirani R. (2017). The elements of statistical learning.
New York: Springer.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{alfa.nb}
}
}

\examples{
x <- as.matrix(iris[, 1:4])
x <- x / rowSums(x)
mod <- alfanb.tune(x, ina = iris[, 5], a = c(0, 0.1, 0.2) )
}

