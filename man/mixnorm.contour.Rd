\name{Contour plot of the Gaussian mixture model in S^2}
\alias{mixnorm.contour}
\title{
Contour plot of the Gaussian mixture model in \eqn{S^2}
}

\description{
Contour plot of the Gaussian mixture model in \eqn{S^2}.
}

\usage{
mixnorm.contour(x, mod)
}

\arguments{
\item{x}{
A matrix with the compositional data.
}
\item{mod}{
An object containing the output of a \code{\link{mix.compnorm}} model.
}
}

\details{
The contour plot of a Gaussian mixture model is plotted. For this you need the data and the fitted model.
}
\value{
A ternary plot with the data and the contour lines of the fitted Gaussian mixture model.
}

\references{
Ryan P. Browne, Aisha ElSherbiny and Paul D. McNicholas (2015). R package mixture: Mixture Models for Clustering and Classification

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
\code{\link{mix.compnorm}, \link{bic.mixcompnorm}, \link{diri.contour}
}
}

\examples{
\dontrun{
x <- as.matrix(iris[, 1:3])
x <- x / rowSums(x)
mod <- mix.compnorm(x, 3, model = "EII")
mixnorm.contour(x, mod)
}
}

\keyword{ contour plot }
\keyword{ Gaussianmixture model }
