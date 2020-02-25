\name{Divergence matrix of compositional data}
\alias{divergence}
\title{
Divergence matrix of compositional data
}

\description{
Divergence matrix of compositional data.
}

\usage{
divergence(x, type = "kullback_leibler", vector = FALSE)
}

\arguments{
\item{x}{
A matrix with the compositional data.
}
\item{type}{
This is either "kullback_leibler" (Kullback-Leibler, which computes the symmetric Kullback-Leibler divergence) or "jensen_shannon" (Jensen-Shannon) divergence.
}
\item{vector}{
For return a vector instead a matrix.
}
}

\details{
The function produces the distance matrix either using the Kullback-Leibler (distance) or the Jensen-Shannon (metric) divergence. The Kullback-Leibler refers to the symmetric Kullback-Leibler divergence.
}

\value{
if the vector argument is FALSE a symmetric matrix with the divergences, otherwise a vector with the divergences.
}

\references{
Endres, D. M. and Schindelin, J. E. (2003). A new metric for probability distributions.
Information Theory, IEEE Transactions on 49, 1858-1860.

Osterreicher, F. and Vajda, I. (2003). A new class of metric divergences on probability
spaces and its applicability in statistics. Annals of the Institute of Statistical
Mathematics 55, 639-653.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{comp.knn}, \link{js.compreg}
}
}

\examples{
x <- as.matrix(iris[1:20, 1:4])
x <- x / rowSums(x)
divergence(x)
}

