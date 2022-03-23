\name{Kernel regression with a numerical response vector or matrix}
\alias{kern.reg}

\title{
Kernel regression with a numerical response vector or matrix
}

\description{
Kernel regression (Nadaraya-Watson estimator) with a numerical response vector or matrix.
}

\usage{
kern.reg(xnew, y, x, h = seq(0.1, 1, length = 10), type = "gauss" )
}

\arguments{
\item{xnew}{
A matrix with the new predictor variables whose compositions are to be predicted.
}
\item{y}{
A numerical vector or a matrix with the response value.
}
\item{x}{
A matrix with the available predictor variables.
}
\item{h}{
The bandwidth value(s) to consider.
}
\item{type}{
The type of kernel to use, "gauss" or "laplace".
}
}

\details{
The Nadaraya-Watson estimator regression is applied.
}

\value{
The fitted values. If a single bandwidth is considered then this is a vector or a matrix, depeding on the nature of the response. If multiple bandwidth values are considered then this is a matrix, if the response is a vector, or a list, if the response is a matrix.
}

\references{
Wand M. P. and Jones M. C. (1994). Kernel smoothing. CRC press.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
  %%  ~~further notes~~
    %}

\seealso{
\code{\link{kernreg.tune}, \link{ice.kernreg}, \link{akern.reg}, \link{aknn.reg}}
}

\examples{
y <- iris[, 1]
x <- iris[, 2:4]
est <- kern.reg(x, y, x, h = c(0.1, 0.2) )
}


