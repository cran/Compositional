\name{The k-NN algorithm for compositional data}
\alias{comp.knn}
\alias{alfa.knn}
\alias{ait.knn}

\title{
The k-NN algorithm for compositional data
}

\description{
The k-NN algorithm for compositional data with and without using the
power transformation.
}

\usage{
comp.knn(xnew, x, ina, a = 1, k = 5,  apostasi = "ESOV", mesos = TRUE)

alfa.knn(xnew, x, ina, a = 1, k = 5, mesos = TRUE,
apostasi = "euclidean", rann = FALSE)

ait.knn(xnew, x, ina, a = 1, k = 5, mesos = TRUE,
apostasi = "euclidean", rann = FALSE)
}

\arguments{
\item{xnew}{
A matrix with the new compositional data whose group is to be predicted. Zeros
are allowed, but you must be carefull to choose strictly positive values
of \eqn{\alpha} or not to set apostasi= "Ait".
}
\item{x}{
A matrix with the available compositional data. Zeros are allowed, but you
must be carefull to choose strictly positive values of \eqn{\alpha} or not
to set apostasi= "Ait".
}
\item{ina}{
A group indicator variable for the available data.
}
\item{a}{
The value of \eqn{\alpha}. As zero values in the compositional data are allowed,
you must be careful to choose strictly positive vcalues of \eqn{\alpha}.
You have the option to put a = NULL. In this case, the xnew and x are
assumed to be the already \eqn{\alpha}-transformed data.
}
\item{k}{
The number of nearest neighbours to consider. It can be a single number or a vector.
}
\item{apostasi}{
The type of distance to use. For the compk.knn this can be one of the following:
"ESOV", "taxicab", "Ait", "Hellinger", "angular" or "CS". See the references for
them. For the alfa.knn this can be either "euclidean" or "manhattan".
}
\item{mesos}{
This is used in the non standard algorithm. If TRUE, the arithmetic mean of the
distances is calulated, otherwise the harmonic mean is used (see details).
}
\item{rann}{
If you have large scale datasets and want a faster k-NN search, you can use
kd-trees implemented in the R package "Rnanoflann". In this case you must set this
argument equal to TRUE. Note however, that in this case, the only available
distance is by default "euclidean".
}
}

\details{
The k-NN algorithm is applied for the compositional data. There are many metrics
and possibilities to choose from. The algorithm finds the k nearest observations
to a new observation and allocates it to the class which appears most times in
the neighbours. It then computes the arithmetic or the harmonic mean of the
distances. The new point is allocated to the class with the minimum distance.
}

\value{
A vector with the estimated groups.
}

\references{
Tsagris, Michail (2014). The k-NN algorithm for compositional data: a revised
approach with and without zero values present. Journal of Data Science, 12(3): 519--534.

Friedman Jerome, Trevor Hastie and Robert Tibshirani (2009). The elements of
statistical learning, 2nd edition. Springer, Berlin

Tsagris Michail, Simon Preston and Andrew T.A. Wood (2016).
Improved classification for compositional data using the
\eqn{\alpha}-transformation. Journal of Classification 33(2): 243--261.

Connie Stewart (2017). An approach to measure distance between compositional
diet estimates containing essential zeros. Journal of Applied Statistics 44(7): 1137--1152.

Clarotto L., Allard D. and Menafoglio A. (2022). A new class of
\eqn{\alpha}-transformations for the spatial analysis of Compositional Data.
Spatial Statistics, 47.

Endres, D. M. and Schindelin, J. E. (2003). A new metric for probability
distributions. Information Theory, IEEE Transactions on 49, 1858--1860.

Osterreicher, F. and Vajda, I. (2003). A new class of metric divergences on
probability spaces and its applicability in statistics.
Annals of the Institute of Statistical Mathematics 55, 639--653.
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
\code{\link{compknn.tune}, \link{alfa.rda}, \link{comp.nb}, \link{alfa.nb}, \link{alfa},
\link{esov}, \link{mix.compnorm}
}
}

\examples{
x <- as.matrix( iris[, 1:4] )
x <- x/ rowSums(x)
ina <- iris[, 5]
mod <- comp.knn(x, x, ina, a = 1, k = 5)
table(ina, mod)
mod2 <- alfa.knn(x, x, ina, a = 1, k = 5)
table(ina, mod2)

}

