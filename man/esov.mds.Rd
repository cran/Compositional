\name{Principal coordinate analysis using the Jensen-Shannon divergence}
\alias{esov.mds}
\title{
Principal coordinate analysis using the Jensen-Shannon divergence
}

\description{
Principal coordinate analysis using the Jensen-Shannon divergence.
}

\usage{
esov.mds(x, k = 2, eig = TRUE)
}

\arguments{
\item{x}{
A matrix with the compositional data. Zero values are allowed.
}
\item{k}{
The maximum dimension of the space which the data are to be represented in. This can be a number between
1 and \eqn{D-1}, where \eqn{D} denotes the number of dimensions.
}
\item{eig}{
Should eigenvalues be returned? The default value is TRUE.
}
}

\details{
The function computes the Jensen-Shannon divergence matrix and then plugs it into the classical
multidimensional scaling function in the "cmdscale" function.
}

\value{
A list with the results of "cmdscale" function.
}

\references{
Aitchison J. (1986). The statistical analysis of compositional data. Chapman & Hall.

Cox, T. F. and Cox, M. A. A. (2001). Multidimensional Scaling. Second edition. Chapman and Hall.

Mardia, K. V., Kent, J. T. and Bibby, J. M. (1979). Chapter 14 of Multivariate Analysis, London: Academic Press.

Tsagris, Michail (2015). A novel, divergence based, regression for compositional data.
Proceedings of the 28th Panhellenic Statistics Conference, 15-18/4/2015, Athens, Greece.
https://arxiv.org/pdf/1511.07600.pdf
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{alfa.mds}, \link{alfa.pca},
}
}

\examples{
x <- as.matrix(iris[, 1:4])
x <- x/ rowSums(x)
a <- esov.mds(x)
}


