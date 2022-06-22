\name{Mixture model selection with the alpha-transformation using BIC}
\alias{bic.alfamixnorm}
\title{
Mixture model selection with the \eqn{\alpha}-transformation using BIC
}

\description{
Mixture model selection with the \eqn{\alpha}-transformation using BIC.
}

\usage{
bic.alfamixnorm(x, G, a = seq(-1, 1, by = 0.1), veo = FALSE, graph = TRUE)
}

\arguments{
\item{x}{
A matrix with compositional data.
}
\item{G}{
A numeric vector with the number of components, clusters, to be considered, e.g. 1:3.
}
\item{a}{
A vector with a grid of values of the power transformation, it has to be between -1 and 1. If zero values are present
it has to be greater than 0. If \eqn{\alpha=0} the isometric log-ratio transformation is applied.
}
\item{veo}{
Stands for "Variables exceed observations". If TRUE then if the number variablesin the model exceeds the number of
observations, but the model is still fitted.
}
\item{graph}{
A boolean variable, TRUE or FALSE specifying whether a graph should be drawn or not.
}
}

\details{
The \eqn{\alpha}-transformation is applied to the compositional data first and then mixtures of multivariate Gaussian
distributions are fitted. BIC is used to decide on the optimal model and number of components.
}

\value{
A list including:

\item{abic}{
A list that contains the matrices of all BIC values for all values of \eqn{\alpha}.
}
\item{optalpha}{
The value of \eqn{\alpha} that leads to the highest BIC.
}
\item{optG}{
The number of components with the highest BIC.
}
\item{optmodel}{
The type of model corresponding to the highest BIC.
}

If graph is set equal to TRUE a plot with the BIC of the best model for each number of components versus the number of components and a list with the results of the Gaussian mixture model for each value of \eqn{\alpha}.
}

\references{
Ryan P. Browne, Aisha ElSherbiny and Paul D. McNicholas (2018). mixture: Mixture Models for Clustering and Classification.
R package version 1.5.

Ryan P. Browne and Paul D. McNicholas (2014). Estimating Common Principal Components in High Dimensions. Advances
in Data Analysis and Classification, 8(2), 217-226.

Aitchison J. (1986). The statistical analysis of compositional data. Chapman & Hall.

Tsagris M.T., Preston S. and Wood A.T.A. (2011). A data-based power transformation for compositional data.
In Proceedings of the 4th Compositional Data Analysis Workshop, Girona, Spain.
https://arxiv.org/pdf/1106.1451.pdf
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{alfa.mix.norm}, \link{mix.compnorm}, \link{mix.compnorm.contour}, \link{rmixcomp}, \link{alfa}, \link{alfa.knn},
\link{alfa.rda}, \link{comp.nb}
}
}

\examples{
\donttest{
x <- as.matrix( iris[, 1:4] )
x <- x/ rowSums(x)
bic.alfamixnorm(x, 1:3, a = c(0.4, 0.5, 0.6), graph = FALSE)
}
}

