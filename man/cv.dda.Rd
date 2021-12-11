\name{Cross-validation for the Dirichlet discriminant analysis}
\alias{cv.dda}
\title{
Cross-validation for the Dirichlet discriminant analysis
}

\description{
Cross-validation for the Dirichlet discriminant analysis.
}

\usage{
cv.dda(x, ina, nfolds = 10, folds = NULL, stratified = TRUE, seed = FALSE)
}
\arguments{
\item{x}{
A matrix with the available data, the predictor variables.
}
\item{ina}{
A vector of data. The response variable, which is categorical (factor is acceptable).
}
\item{folds}{
A list with the indices of the folds.
}
\item{nfolds}{
The number of folds to be used. This is taken into consideration only if "folds" is NULL.
}
\item{stratified}{
Do you want the folds to be selected using stratified random sampling? This preserves the analogy of the samples of each group. Make this TRUE if you wish.
}
\item{seed}{
If you set this to TRUE, the same folds will be created every time.
}
}

\details{
This function estimates the performance of the Dirichlet discriminant analysis via k-fold cross-validation.
}

\value{
A list including:
\item{percent}{
The percentage of correct classification
}
\item{runtime}{
The duration of the cross-validation proecdure.
}
}

\references{
Friedman J., Hastie T. and Tibshirani R. (2017). The elements of statistical learning.
New York: Springer.

Thomas P. Minka (2003). Estimating a Dirichlet distribution.
http://research.microsoft.com/en-us/um/people/minka/papers/dirichlet/minka-dirichlet.pdf
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{dda}, \link{alfanb.tune}, \link{alfarda.tune}, \link{compknn.tune}, \link{cv.compnb}
}
}

\examples{
x <- as.matrix(iris[, 1:4])
x <- x / rowSums(x)
mod <- cv.dda(x, ina = iris[, 5] )
}

