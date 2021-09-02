\name{Mixture model selection via BIC}
\alias{bic.mixcompnorm}
\title{
Mixture model selection via BIC
}

\description{
Mixture model selection via BIC.
}

\usage{
bic.mixcompnorm(x, G, type = "alr", veo = FALSE, graph = TRUE)
}

\arguments{
\item{x}{
A matrix with compositional data.
}
\item{G}{
A numeric vector with the number of components, clusters, to be considered.
}
\item{type}{
The type of trasformation to be used, either the additive log-ratio ("alr"), the isometric log-ratio ("ilr") or the pivot coordinate ("pivot") transformation.
}
\item{veo}{
Stands for "Variables exceed observations". If TRUE then if the number variablesin the model exceeds the number of observations, but the model is still fitted.
}
\item{graph}{
A boolean variable, TRUE or FALSE specifying whether a graph should be drawn or not.
}
}

\details{
The alr or the ilr-transformation is applied to the compositional data first and then mixtures of multivariate Gaussian distributions are fitted. BIC is used to decide on the optimal model and number of components.
}

\value{
A plot with the BIC of the best model for each number of components versus the number of components.
A list including:
\item{mod}{
A message informing the user about the best model.
}
\item{BIC}{
The BIC values for every possible model and number of components.
}
}

\references{
Ryan P. Browne, Aisha ElSherbiny and Paul D. McNicholas (2018). mixture: Mixture Models for Clustering and Classification. R package version 1.5.

Ryan P. Browne and Paul D. McNicholas (2014). Estimating Common Principal Components in High Dimensions. Advances in Data Analysis and Classification, 8(2), 217-226.

Aitchison J. (1986). The statistical analysis of compositional data. Chapman & Hall.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{mix.compnorm}, \link{mixnorm.contour}, \link{rmixcomp}, \link{bic.alfamixnorm}
}
}

\examples{
\dontrun{
x <- as.matrix( iris[, 1:4] )
x <- x/ rowSums(x)
bic.mixcompnorm(x, 1:3, type = "alr", graph = FALSE)
bic.mixcompnorm(x, 1:3, type = "ilr", graph = FALSE)
}
}

\keyword{ mixtures of normal distributions }
\keyword{ model selection }
