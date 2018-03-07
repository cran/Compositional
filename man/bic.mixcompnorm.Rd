\name{Mixture model selection via BIC}
\alias{bic.mixcompnorm}
\title{
Mixture model selection via BIC
}

\description{
Mixture model selection via BIC.
}

\usage{
bic.mixcompnorm(x, A, type = "alr")
}

\arguments{
\item{x}{
A matrix with compositional data.
}
\item{A}{
The maximum number of components, clusters, to be considered.
}
\item{type}{
The type of trasformation to be used, either additive log-ratio ("alr") or the isometric log-ratio ("ilr").
}
}

\details{
The alr or the ilr-transformation is applied to the compositional data first and then mixtures of multivariate Gaussian distributions are fitted.
BIC is used to decide on the optimal model and number of components.
}

\value{
a plot with the BIC of the best model for each number of components versus the number of components.
A list including:
\item{mod}{
A message informing the user about the best model.
}
\item{BIC}{
The BIC values for every possible model and number of components.
}
}

\references{
Ryan P. Browne, Aisha ElSherbiny and Paul D. McNicholas (2015). mixture: Mixture Models for Clustering and Classification. R package version 1.4.

Ryan P. Browne and Paul D. McNicholas (2014). Estimating Common Principal Components in High Dimensions. Advances in Data Analysis and Classification, 8(2), 217-226.

Aitchison J. (1986). The statistical analysis of compositional data. Chapman \& Hall.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <gioathineou@gmail.com>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{mix.compnorm}, \link{mixnorm.contour}, \link{rmixcomp}
}
}

\examples{
library(MASS)
x <- as.matrix( iris[, 1:4] )
x <- x/ rowSums(x)
bic.mixcompnorm(x, 6, type = "alr")
bic.mixcompnorm(x, 6, type = "ilr")
}

\keyword{ mixtures of normal distributions }
\keyword{ model selection }
