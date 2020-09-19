\name{The alpha-transformation}
\alias{alfa}
\alias{alef}
\title{
The \eqn{\alpha}-transformation
}

\description{
The \eqn{\alpha}-transformation.
}

\usage{
alfa(x, a, h = TRUE)
alef(x, a)
}

\arguments{
\item{x}{
A matrix with the compositional data.
}
\item{a}{
The value of the power transformation, it has to be between -1 and 1. If zero values are present it has to
be greater than 0. If \eqn{\alpha=0} the isometric log-ratio transformation is applied.
}
\item{h}{
A boolean variable. If is TRUE (default value) the multiplication with the Helmert sub-matrix will take place.
When \eqn{\alpha=0} and h = FALSE, the result is the centred log-ratio transformation (Aitchison, 1986).
In general, when h = FALSE the resulting transformation maps the data onto a singualr space. The sum of the vectors is equal to 0.
Hence, from the simplex constraint the data go to another constraint.
}
}

\details{
The \eqn{\alpha}-transformation is applied to the compositional data. The command "alef" is the same as
"alfa(x, a, h = FALSE)", but reurns a different element as well and is necessary for the functions \code{\link{a.est}}, \code{\link{a.mle}} and \code{\link{alpha.mle}}.
}

\value{
A list including:
\item{sa}{
The logarithm of the Jacobian determinant of the \eqn{\alpha}-transformation. This is used in the "profile"
function to speed up the computations.
}
\item{sk}{
If the "alef" was called, this will return the sum of the \eqn{\alpha}-power transformed data, prior to
being normalised to sum to 1. If \eqn{\alpha=0}, this will not be returned.
}
\item{aff}{
The \eqn{\alpha}-transformed data.
}
}

\references{
Tsagris M.T., Preston S. and Wood A.T.A. (2011). A data-based power transformation for compositional data. In Proceedings of the 4th Compositional Data Analysis Workshop, Girona, Spain.
https://arxiv.org/pdf/1106.1451.pdf

Tsagris Michail and Stewart Connie, (2020). A folded model for compositional data analysis. Australian and New Zealnd Journal of Statistics (to appear).
https://arxiv.org/pdf/1802.07330.pdf

Aitchison J. (1986). The statistical analysis of compositional data. Chapman & Hall.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} and Giorgos Athineou <gioathineou@gmail.com>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{alfainv}, \link{alfa.profile}, \link{alfa.tune}
\link{a.est}, \link{alpha.mle}
}
}

\examples{
library(MASS)
x <- as.matrix(fgl[, 2:9])
x <- x / rowSums(x)
y1 <- alfa(x, 0.2)$aff
y2 <- alfa(x, 1)$aff
rbind( colMeans(y1), colMeans(y2) )
y3 <- alfa(x, 0.2)$aff
dim(y1)  ;  dim(y3)
rowSums(y1)
rowSums(y3)
}

