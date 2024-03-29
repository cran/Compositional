\name{Inverse of the alpha-transformation}
\alias{alfainv}
\title{
Inverse of the \eqn{\alpha}-transformation
}

\description{
The inverse of the \eqn{\alpha}-transformation.
}

\usage{
alfainv(x, a, h = TRUE)
}

\arguments{
  \item{x}{
  A matrix with Euclidean data. However, they must lie within the feasible,
  acceptable space. See references for more information.
}
  \item{a}{
  The value of the power transformation, it has to be between -1 and 1.
  If zero values are present it has to be greater than 0. If \eqn{\alpha=0},
  the inverse of the isometric log-ratio transformation is applied.
}
\item{h}{
  If h = TRUE this means that the multiplication with the Helmer sub-matrix
  will take place. It is set to TRUe by default.
}
}

\details{
The inverse of the \eqn{\alpha}-transformation is applied to the data.
If the data lie outside the \eqn{\alpha}-space, NAs will be returned for
some values.
}

\value{
A matrix with the pairwise distances.
}

\references{
Tsagris M. and Stewart C. (2022). A Review of Flexible Transformations for Modeling Compositional Data. 
In Advances and Innovations in Statistics and Data Science, pp. 225--234.
https://link.springer.com/chapter/10.1007/978-3-031-08329-7_10

Tsagris M.T., Preston S. and Wood A.T.A. (2016). Improved classification for
compositional data using the \eqn{\alpha}-transformation. 
Journal of Classification 33(2): 243--261.
https://arxiv.org/pdf/1506.04976v2.pdf

Tsagris M.T., Preston S. and Wood A.T.A. (2011). A data-based power
transformation for compositional data.
In Proceedings of the 4th Compositional Data Analysis Workshop, Girona, Spain.
https://arxiv.org/pdf/1106.1451.pdf
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
\code{\link{alfa}, \link{alfadist}
}
}

\examples{
library(MASS)
x <- as.matrix(fgl[1:10, 2:9])
x <- x / rowSums(x)
y <- alfa(x, 0.5)$aff
alfainv(y, 0.5)
}

\keyword{ inverse transformation }
