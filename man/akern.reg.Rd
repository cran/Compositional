\name{The alpha-kernel regression for compositional response data}
\alias{akern.reg}

\title{
The \eqn{\alpha}-kernel regression for compositional response data
}

\description{
The \eqn{\alpha}-kernel regression for compositional response data.
}

\usage{
akern.reg( xnew, y, x, a = seq(0.1, 1, by = 0.1), h = seq(0.1, 1, length = 10) )
}

\arguments{
\item{xnew}{
A matrix with the new predictor variables whose compositions are to be predicted.
}
\item{y}{
A matrix with the compositional response data. Zeros are allowed.
}
\item{x}{
A matrix with the available predictor variables.
}
\item{a}{
The value(s) of \eqn{\alpha}. Either a single value or a vector of values. As zero values in the compositional data are allowed, you must be careful to choose strictly positive vcalues of \eqn{\alpha}. However, if negative values are passed, the positive ones are used only.
}
\item{h}{
The bandwidth value(s) to consider.
}
}

\details{
The \eqn{\alpha}-kernel regression for compositional response variables is applied.
}

\value{
A list with the estimated compositional response data for each value of \eqn{\alpha} and h.
}

\references{
Michail Tsagris, Abdulaziz Alenazi and Connie Stewart (2021).
Non-parametric regression models for compositional data.
https://arxiv.org/pdf/2002.05137.pdf
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
  %%  ~~further notes~~
    %}

\seealso{
\code{\link{akernreg.tune}, \link{aknn.reg}, \link{aknnreg.tune}, \link{alfa.reg}, \link{comp.ppr}, \link{comp.reg}, \link{kl.compreg}
}
}

\examples{
y <- as.matrix( iris[, 1:3] )
y <- y / rowSums(y)
x <- iris[, 4]
mod <- akern.reg( x, y, x, a = c(0.4, 0.5), h = c(0.1, 0.2) )
}


