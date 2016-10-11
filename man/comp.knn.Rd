\name{The k-NN algorithm for compositional data}
\alias{comp.knn}
\alias{alfa.knn}

\title{
The k-NN algorithm for compositional data
}

\description{
The k-NN algorithm for compositional data with and without using the power transformation.
}

\usage{
comp.knn(xnew, x, ina, a = 1, k = 5, type = "S", apostasi = "ESOV", mesos = TRUE)

alfa.knn(xnew, x, ina, a = 1, k = 5, type = "S", mesos = TRUE)

}

\arguments{
  \item{xnew}{
  A matrix with the new compositional data whose group is to be predicted. Zeros are allowed, but you must be carefull to choose strictly positive vcalues of \eqn{\alpha} or not to set apostasi= "Ait".
}
  \item{x}{
  A matrix with the available compositional data. Zeros are allowed, but you must be carefull to choose strictly positive vcalues of \eqn{\alpha} or not to set apostasi= "Ait".
}
  \item{ina}{
  A group indicator variable for the avaiable data.
}
  \item{a}{
The value of \eqn{\alpha}. As zero values in the compositional data are allowed, you must be carefull to choose strictly positive vcalues of \eqn{\alpha}.
}
  \item{k}{
 The number of nearest neighbours to consider.
}
  \item{type}{
  This can be either "S" for the standard k-NN or "NS" for the non standard (see details).
  }
  \item{apostasi}{
  The type of distance to use. "ESOV", "taxicab", "Ait", "Hellinger", "angular" or "CS". See the references for them.
}
  \item{mesos}{
 This is used in the non standard algorithm. If TRUE, the arithmetic mean of the distances is calulated, otherwise the harmonic mean is used (see details).
}
}

\details{
The k-NN algorithm is applied for the compositional data. There are many metrics and possibilities to choose from. The standard algorithm finds the k nearest observations to a new observation and allocates it to the class which appears most times in the neighbours. The non standard algorithm is slower but perhaps more accurate. For every group is finds the k nearest neighbours to the new observation. It then computes the arithmetic or the harmonic mean of the distances. The new point is allocated to the class with the minimum distance.
}

\value{
A vector with the estimated groups.
}

\references{
Tsagris, Michail (2014). The k-NN algorithm for compositional data: a revised approach with and without zero values present. Journal of Data Science, 12(3): 519-534.

Friedman Jerome, Trevor Hastie and Robert Tibshirani (2009). The elements of statistical learning, 2nd edition. Springer, Berlin

Tsagris Michail, Simon Preston and Andrew TA Wood (2016). Improved classification for compositional data using the \eqn{\alpha}-transformation. Journal of classification (to appear).
http://arxiv.org/pdf/1106.1451.pdf

Connie Stewart (2016). An approach to measure distance between compositional diet estimates containing essential zeros.
Journal of Applied Statistics. http://www.tandfonline.com/doi/full/10.1080/02664763.2016.1193846
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{compknn.tune}, \link{rda}, \link{alfa}
}
}

\examples{
x <- as.matrix( iris[, 1:4] )
ina <- iris[, 5]
mod <- comp.knn(x, x, ina, a = 1, k = 5)
table(ina, mod)
mod2 <- alfa.knn(x, x, ina, a = 1, k = 5)
table(ina, mod2)

}

\keyword{ Regularised discriminant analysis }

