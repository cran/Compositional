\name{Tuning of the k-NN algorithm for compositional data}
\alias{compknn.tune}
\alias{alfaknn.tune}
\alias{aitknn.tune}

\title{
Tuning of the k-NN algorithm for compositional data
}

\description{
Tuning of the k-NN algorithm for compositional data with and without using the
power or the \eqn{\alpha}-transformation. In addition, estimation of the rate
of correct classification via K-fold cross-validation.
}

\usage{
compknn.tune(x, ina, nfolds = 10, k = 2:5, mesos = TRUE,
a = seq(-1, 1, by = 0.1), apostasi = "ESOV", folds = NULL,
stratified = TRUE, seed = NULL, graph = FALSE)

alfaknn.tune(x, ina, nfolds = 10, k = 2:5, mesos = TRUE,
a = seq(-1, 1, by = 0.1), apostasi = "euclidean", rann = FALSE,
folds = NULL, stratified = TRUE, seed = NULL, graph = FALSE)

aitknn.tune(x, ina, nfolds = 10, k = 2:5, mesos = TRUE,
a = seq(-1, 1, by = 0.1), apostasi = "euclidean", rann = FALSE,
folds = NULL, stratified = TRUE, seed = NULL, graph = FALSE)
}

\arguments{
\item{x}{
A matrix with the available compositional data. Zeros are allowed, but you
must be careful to choose strictly positive values of \eqn{\alpha} or not
to set apostasi= "Ait".
}
\item{ina}{
A group indicator variable for the available data.
}
\item{nfolds}{
The number of folds to be used. This is taken into consideration only if
the folds argument is not supplied.
}
\item{k}{
A vector with the nearest neighbours to consider.
}
\item{mesos}{
This is used in the non standard algorithm. If TRUE, the arithmetic mean of
the distances is calculated, otherwise the harmonic mean is used (see details).
}
\item{a}{
A grid of values of \eqn{\alpha} to be used only if the distance chosen allows
for it.
}
\item{apostasi}{
The type of distance to use. For the compk.knn this can be one of the following:
"ESOV", "taxicab", "Ait", "Hellinger", "angular" or "CS". See the references
for them. For the alfa.knn this can be either "euclidean" or "manhattan".
}
\item{rann}{
If you have large scale datasets and want a faster k-NN search, you can use
kd-trees implemented in the R package "RANN". In this case you must set this
argument equal to TRUE. Note however, that in this case, the only available
distance is by default "euclidean".
}
\item{folds}{
If you have the list with the folds supply it here. You can also leave it NULL
and it will create folds.
}
\item{stratified}{
Do you want the folds to be created in a stratified way? TRUE or FALSE.
}
\item{seed}{
You can specify your own seed number here or leave it NULL.
}
\item{graph}{
If set to TRUE a graph with the results will appear.
}
}

\details{
The k-NN algorithm is applied for the compositional data. There are many metrics
and possibilities to choose from. The algorithm finds the k nearest
observations to a new observation and allocates it to the class which appears
most times in the neighbours.
}

\value{
A list including:
\item{per}{
A matrix or a vector (depending on the distance chosen) with the averaged over
all folds rates of correct classification for all hyper-parameters
(\eqn{\alpha} and k).
}
\item{performance}{
The estimated rate of correct classification.
}
\item{best_a}{
The best value of \eqn{\alpha}. This is returned for "ESOV" and "taxicab" only.
}
\item{best_k}{
The best number of nearest neighbours.
}
\item{runtime}{
The run time of the cross-validation procedure.
}
}

\references{
Tsagris, Michail (2014). The k-NN algorithm for compositional data:
a revised approach with and without zero values present.
Journal of Data Science, 12(3): 519-534.
https://arxiv.org/pdf/1506.05216.pdf

Friedman Jerome, Trevor Hastie and Robert Tibshirani (2009).
The elements of statistical learning,
2nd edition. Springer, Berlin

Tsagris M., Preston S. and Wood A.T.A. (2016). Improved classification for
compositional data using the \eqn{\alpha}-transformation.
Journal of classification, 33(2):243-261.
http://arxiv.org/pdf/1106.1451.pdf

Connie Stewart (2017). An approach to measure distance between
compositional diet estimates containing essential zeros.
Journal of Applied Statistics 44(7): 1137-1152.

Clarotto L., Allard D. and Menafoglio A. (2022).
A new class of \eqn{\alpha}-transformations for the spatial analysis
of Compositional Data. Spatial Statistics, 47.

Endres, D. M. and Schindelin, J. E. (2003). A new metric for probability distributions.
Information Theory, IEEE Transactions on 49, 1858-1860.

Osterreicher, F. and Vajda, I. (2003). A new class of metric divergences on
probability spaces and its applicability in statistics.
Annals of the Institute of Statistical Mathematics 55, 639-653.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}
and Giorgos Athineou <gioathineou@gmail.com>.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{comp.knn}, \link{alfarda.tune}, \link{cv.dda}, \link{cv.compnb}
}
}

\examples{
x <- as.matrix(iris[, 1:4])
x <- x/ rowSums(x)
ina <- iris[, 5]
mod1 <- compknn.tune(x, ina, a = seq(1, 1, by = 0.1) )
mod2 <- alfaknn.tune(x, ina, a = seq(-1, 1, by = 0.1) )
}

\keyword{ Regularised discriminant analysis }

