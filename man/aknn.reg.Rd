\name{The alpha-k-NN regression for compositional response data}
\alias{aknn.reg}

\title{
The \eqn{\alpha}-k-NN regression for compositional response data
}

\description{
The \eqn{\alpha}-k-NN regression for compositional response data.
}

\usage{
aknn.reg(xnew, y, x, a = seq(0.1, 1, by = 0.1), k = 2:10,
apostasi = "euclidean", rann = FALSE)
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
The value(s) of \eqn{\alpha}. Either a single value or a vector of values.
As zero values in the compositional data are allowed, you must be careful
to choose strictly positive vcalues of \eqn{\alpha}. However, if negative
values are passed, the positive ones are used only.
}
\item{k}{
The number of nearest neighbours to consider. It can be a single number or a vector.
}
\item{apostasi}{
The type of distance to use, either "euclidean" or "manhattan".
}
\item{rann}{
If you have large scale datasets and want a faster k-NN search, you can use kd-trees implemented in the R package "RANN". In this case you must set this argument equal to TRUE. Note however, that in this case, the only available distance is by default "euclidean".
}
}

\details{
The \eqn{\alpha}-k-NN regression for compositional response variables is applied.
}

\value{
A list with the estimated compositional response data for each value of \eqn{\alpha} and k.
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
\code{\link{aknnreg.tune}, \link{akern.reg}, \link{alfa.reg}, \link{comp.ppr}, \link{comp.reg}, \link{kl.compreg}
}
}

\examples{
y <- as.matrix( iris[, 1:3] )
y <- y / rowSums(y)
x <- iris[, 4]
mod <- aknn.reg(x, y, x, a = c(0.4, 0.5), k = 2:3, apostasi = "euclidean")
}


