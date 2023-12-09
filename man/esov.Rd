\name{The ESOV-distance}
\alias{esov}
\alias{esova}
\alias{es}
\title{
The ESOV-distance
}

\description{
The ESOV-distance.
}

\usage{
esov(x)
esova(xnew, x)
es(x1, x2)
}

\arguments{
\item{x}{
A matrix with compositional data.
}
\item{xnew}{
A matrix or a vector with new compositional data.
}
\item{x1}{
A vector with compositional data.
}

\item{x2}{
A vector with compositional data.
}
}

\details{
The ESOV distance is calculated.
}

\value{
For "esov()" a matrix including the pairwise distances of all observations or the distances between xnew and x.

For "esova()" a matrix including the pairwise distances of all observations or the distances between xnew and x.

For "es()" a number, the ESOV distance between x1 and x2.
}


\references{
Tsagris, Michail (2014). The k-NN algorithm for compositional data: a revised approach with and without zero values present. Journal of Data Science, 12(3): 519-534.

Endres, D. M. and Schindelin, J. E. (2003). A new metric for probability distributions. Information Theory, IEEE Transactions on 49, 1858-1860.

Osterreicher, F. and Vajda, I. (2003). A new class of metric divergences on probability spaces and its applicability in statistics. Annals of the Institute of Statistical Mathematics 55, 639-653.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{alfadist}, \link{comp.knn}, \link{js.compreg}
}
}

\examples{
library(MASS)
x <- as.matrix(fgl[1:20, 2:9])
x <- x / rowSums(x)
esov(x)
}

