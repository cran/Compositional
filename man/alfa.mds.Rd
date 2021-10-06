\name{Principal coordinate analysis using the alpha-distance}
\alias{alfa.mds}
\title{
Principal coordinate analysis using the \eqn{\alpha}-distance
}

\description{
Principal coordinate analysis using the \eqn{\alpha}-distance.
}

\usage{
alfa.mds(x, a, k = 2, eig = TRUE)
}

\arguments{
\item{x}{
A matrix with the compositional data. Zero values are allowed.
}
\item{a}{
The value of a. In case of zero values in the data it has to be greater than 1.
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
The function computes the \eqn{\alpha}-distance matrix and then plugs it into the classical
multidimensional scaling function in the "cmdscale" function.
}

\value{
A list with the results of "cmdscale" function.
}

\references{
Aitchison J. (1986). The statistical analysis of compositional data. Chapman & Hall.

Cox, T. F. and Cox, M. A. A. (2001). Multidimensional Scaling. Second edition. Chapman and Hall.

Mardia, K. V., Kent, J. T. and Bibby, J. M. (1979). Chapter 14 of Multivariate Analysis, London: Academic Press.

Tsagris M.T., Preston S. and Wood A.T.A. (2011). A data-based power transformation for compositional data.
In Proceedings of the 4th Compositional Data Analysis Workshop, Girona, Spain.
https://arxiv.org/pdf/1106.1451.pdf
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{esov.mds}, \link{alfa.pca},
}
}

\examples{
  x <- as.matrix(iris[, 1:4])
  x <- x/ rowSums(x)
  a <- esov.mds(x)
}


