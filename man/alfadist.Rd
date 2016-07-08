\name{The alpha-distance}
\alias{alfadist}
\title{
The \eqn{\alpha}-distance
}

\description{
This is the Euclidean distance after the \eqn{\alpha}-transformation has been applied.
}

\usage{
alfadist(x, a)
}

\arguments{
  \item{x}{
  A matrix with the compositional data.
}
  \item{a}{
  The value of the power transformation, it has to be between -1 and 1. If zero values are present it has to be greater than 0. If \eqn{\alpha=0} the isometric log-ratio transformation is applied and the solution exists in a closed form, since it the classical mutivariate regression.
}
}

\details{
The \eqn{\alpha}-transformation is applied to the compositional data first and then the Euclidean distance is calculated.
}

\value{
A matrix including:
The pairwise distances of all observations.
}


\references{
Tsagris M.T., Preston S. and Wood A.T.A. (2016). Improved classification for compositional data using the \eqn{\alpha}-transformation. Journal of Classification (to appear). http://arxiv.org/pdf/1506.04976v2.pdf

Tsagris M.T., Preston S. and Wood A.T.A. (2011). A data-based power transformation for compositional data. In Proceedings of the 4th Compositional Data Analysis Workshop, Girona, Spain. http://arxiv.org/pdf/1106.1451.pdf
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{alfa}, \link{alfainv}, \link{alfa.reg}
}
}

\examples{
library(MASS)
x <- fgl[1:20, 2:9]
alfadist(x, 0.1)
alfadist(x, 1)
}

\keyword{ \eqn{\alpha}-transformation }
\keyword{ Euclidean distance }
