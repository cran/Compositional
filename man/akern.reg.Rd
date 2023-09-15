\name{The alpha-kernel regression with compositional response data}
\alias{akern.reg}

\title{
The \eqn{\alpha}-kernel regression with compositional response data
}

\description{
The \eqn{\alpha}-kernel regression with compositional response data.
}

\usage{
akern.reg( xnew, y, x, a = seq(0.1, 1, by = 0.1),
h = seq(0.1, 1, length = 10), type = "gauss" )
}

\arguments{
\item{xnew}{
A matrix with the new predictor variables whose compositions are to
be predicted.
}
\item{y}{
A matrix with the compositional response data. Zeros are allowed.
}
\item{x}{
A matrix with the available predictor variables.
}
\item{a}{
The value(s) of \eqn{\alpha}. Either a single value or a vector of values.
As zero values in the compositional data are allowed, you must be careful
to choose strictly positive vcalues of \eqn{\alpha}. However, if negative
values are passed, the positive ones are used only.
}
\item{h}{
The bandwidth value(s) to consider.
}
\item{type}{
The type of kernel to use, "gauss" or "laplace".
}
}

\details{
The \eqn{\alpha}-kernel regression for compositional response variables is
applied.
}

\value{
A list with the estimated compositional response data for each value of
\eqn{\alpha} and h.
}

\references{
Tsagris M., Alenazi A. and Stewart C. (2023).
Flexible non-parametric regression models for compositional response data with zeros.
Statistics and Computing, 33(5): 1--17.
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
\code{\link{akernreg.tune}, \link{aknn.reg}, \link{aknnreg.tune},
\link{alfa.reg}, \link{comp.ppr}, \link{comp.reg}, \link{kl.compreg}
}
}

\examples{
y <- as.matrix( iris[, 1:3] )
y <- y / rowSums(y)
x <- iris[, 4]
mod <- akern.reg( x, y, x, a = c(0.4, 0.5), h = c(0.1, 0.2) )
}


