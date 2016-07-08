\name{Tuning of the k-NN algorithm for compositional data}
\alias{compknn.tune}
\alias{alfaknn.tune}
\title{
Tuning of the he k-NN algorithm for compositional data
}

\description{
Tuning of the k-NN algorithm for compositional data with and without using the power or the \eqn{\alpha}-transformation. In addition, estimation of the rate of correct classification via M-fold cross-validation.
}

\usage{
compknn.tune(x, ina, M = 10, A = 5, type= "S", mesos = TRUE,
a = seq(-1, 1, by = 0.1), apostasi = "ESOV", mat = NULL, graph = FALSE)

alfaknn.tune(x, ina, M = 10, A = 5, type = "S", mesos = TRUE,
a = seq(-1, 1, by = 0.1), mat = NULL, graph = FALSE)
}

\arguments{
\item{x}{
A matrix with the available compositional data. Zeros are allowed, but you must be carefull to choose strictly positive vcalues of \eqn{\alpha} or not to set apostasi= "Ait".
}
\item{ina}{
A group indicator variable for the avaiable data.
}
\item{M}{
The number of folds to be used. This is taken into consideration only if the matrix "mat" is not supplied.
}
\item{A}{
The maximum number of nearest neighbours to consider. Note that the 1 neasrest neighbour is not used.
}
\item{type}{
This can be either "S" for the standard k-NN or "NS" for the non standard (see details).
}
\item{mesos}{
This is used in the non standard algorithm. If TRUE, the arithmetic mean of the distances is calulated, otherwise the harmonic mean is used (see details).
}
\item{a}{
A grid of values of \eqn{\alpha} to be used only if the distance chosen allows for it.
}
\item{apostasi}{
The type of distance to use. "ESOV", "taxicab", "Ait", "Hellinger", "angular" or "CS". See the references for them.
}
\item{mat}{
You can specify your own folds by giving a mat, where each column is a fold. Each column contains indices of the observations.  You can also leave it NULL and it will create folds.
}
\item{graph}{
If set to TRUE a graph with the results will appear.
}
}

\details{
The k-NN algorithm is applied for the compositional data. There are many metrics and possibilities to choose from. The standard algorithm finds the k nearest observations to a new observation and allocates it to the class which appears most times in the neighbours. The non standard algorithm is slower but perhaps more accurate. For every group is finds the k nearest neighbours to the new observation. It then computes the arithmetic or the harmonic mean of the distances. The new point is allocated to the class with the minimum distance.
}

\value{
A list including:
\item{ela}{
A matrix or a vector (depending on the distance chosen) with the averaged over all folds rates of correct classification for all hyper-parameters (\eqn{\alpha} and k).
}
\item{performance}{
The bias corrected estimated rate of correct classification along with the estimated bias.
}
\item{best_a}{
The best value of \eqn{\alpha}. This is returned for "ESOV" and "taxicab" only.
}
\item{best_k}{
The best number of nearest neighbours.
}
\item{runtime}{
The ru time of the cross-validation procedure.
}
}

\references{
Tsagris, Michail (2014). The k-NN algorithm for compositional data: a revised approach with and without zero values present. Journal of Data Science, 12(3): 519-534.

Friedman Jerome, Trevor Hastie and Robert Tibshirani (2009). The elements of statistical learning, 2nd edition. Springer, Berlin

Tsagris Michail, Simon Preston and Andrew TA Wood (2016).  Improved classification for compositional data using the \eqn{\alpha}-transformation. Journal of classification (to appear).
http://arxiv.org/pdf/1106.1451.pdf

Connie Stewart (2016). An approach to measure distance between compositional diet estimates containing essential zeros.
Journal of Applied Statistics.  http://www.tandfonline.com/doi/full/10.1080/02664763.2016.1193846
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{comp.knn}, \link{rda}, \link{alfa}
}
}

\examples{
x <- iris[, 1:4]
ina <- iris[, 5]
mod1 <- compknn.tune(x, ina, a = seq(-0.1, 0.1, by = 0.1) )
mod2 <- alfaknn.tune(x, ina, a = seq(-0.1, 0.1, by = 0.1) )
}

\keyword{ Regularised discriminant analysis }

