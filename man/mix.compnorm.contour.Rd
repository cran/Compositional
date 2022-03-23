\name{Contour plot of the Gaussian mixture model in S^2}
\alias{mix.compnorm.contour}
\title{
Contour plot of the Gaussian mixture model in \eqn{S^2}
}

\description{
Contour plot of the Gaussian mixture model in \eqn{S^2}.
}

\usage{
mix.compnorm.contour(mod, type = "alr", n = 100, x = NULL, cont.line = FALSE)
}

\arguments{
\item{mod}{
An object containing the output of a \code{\link{mix.compnorm}} model.
}
\item{type}{
The type of trasformation used, either the additive log-ratio ("alr"), the isometric log-ratio
("ilr") or the pivot coordinate ("pivot") transformation.
}
\item{n}{
The number of grid points to consider over which the density is calculated.
}
\item{x}{
A matrix with the compositional data.
}
\item{cont.line}{
Do you want the contour lines to appear? If yes, set this TRUE.
}
}

\details{
The contour plot of a Gaussian mixture model is plotted. For this you need the (fitted) model.
}
\value{
A ternary plot with the data and the contour lines of the fitted Gaussian mixture model.
}

\references{
Ryan P. Browne, Aisha ElSherbiny and Paul D. McNicholas (2015). R package mixture: Mixture Models for
Clustering and Classification

Aitchison J. (1986). The statistical analysis of compositional data. Chapman & Hall.
}

\author{
Michail Tsagris and Christos Adam.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} and
Christos Adam \email{pada4m4@gmail.com}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{mix.compnorm}, \link{bic.mixcompnorm}, \link{diri.contour}
}
}

\examples{
x <- as.matrix(iris[, 1:3])
x <- x / rowSums(x)
mod <- mix.compnorm(x, 3, model = "EII")
mix.compnorm.contour(mod, "alr")
}

\keyword{ contour plot }
\keyword{ Gaussianmixture model }
