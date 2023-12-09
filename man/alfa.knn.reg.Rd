\name{The alpha-k-NN regression with compositional predictor variables}
\alias{alfa.knn.reg}

\title{
The \eqn{\alpha}-k-NN regression with compositional predictor variables
}

\description{
The \eqn{\alpha}-k-NN regression with compositional predictor variables.
}

\usage{
alfa.knn.reg(xnew, y, x, a = 1, k = 2:10, apostasi = "euclidean", method = "average")
}

\arguments{
\item{xnew}{
A matrix with the new compositional predictor variables whose response is to be predicted. Zeros are allowed.
}
\item{y}{
The response variable, a numerical vector.
}
\item{x}{
A matrix with the available compositional predictor variables. Zeros are allowed.
}
\item{a}{
A single value of \eqn{\alpha}. As zero values in the compositional data are allowed, you must be careful to choose strictly positive vcalues of \eqn{\alpha}. If negative values are passed, the positive ones are used only. If the data are already alpha-transformed, you can make this NULL.
}
\item{k}{
The number of nearest neighbours to consider. It can be a single number or a vector.
}
\item{apostasi}{
The type of distance to use, either "euclidean" or "manhattan".
}
\item{method}{
If you want to take the average of the reponses of the k closest observations, type "average".
For the median, type "median" and for the harmonic mean, type "harmonic".
}
}

\details{
The \eqn{\alpha}-k-NN regression with compositional predictor variables is applied.
}

\value{
A matrix with the estimated response data for each value of k.
}

\references{
Tsagris M., Alenazi A. and Stewart C. (2023).
Flexible non-parametric regression models for compositional response data with zeros.
Statistics and Computing, 33(106).

https://link.springer.com/article/10.1007/s11222-023-10277-5
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{aknn.reg}, \link{alfa.knn}, \link{alfa.pcr}, \link{alfa.ridge}
}
}

\examples{
library(MASS)
x <- as.matrix(fgl[, 2:9])
x <- x / rowSums(x)
y <- fgl[, 1]
mod <- alfa.knn.reg(x, y, x, a = 0.5, k = 2:4)
}


