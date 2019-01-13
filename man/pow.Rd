\name{Power operation}
\alias{pow}
\title{
Power operation
}

\description{
Power operation.
}

\usage{
pow(x, a)
}

\arguments{
\item{x}{
A matrix with the compositional data.
}
\item{a}{
Either a vector with numbers of a single number.
}
}

\details{
This is the power operation defined by Aitchison (1986). It is also the starting point of the \eqn{\alpha}-transformation.
}

\value{
A matrix with the power transformed compositional data.
}

\references{
Aitchison J. (1986). The statistical analysis of compositional data. Chapman & Hall.

Tsagris M.T., Preston S. and Wood A.T.A. (2011). A data-based power transformation for compositional data.
In Proceedings of the 4th Compositional Data Analysis Workshop, Girona, Spain.
http://arxiv.org/pdf/1106.1451.pdf
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{perturbation}, \link{alfa} }
}

\examples{
x <- as.matrix(iris[1:15, 1:4])
a <- runif(1)
pow(x, a)
}

