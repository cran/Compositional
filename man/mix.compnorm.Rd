\name{Gaussian mixture models for compositional data}
\alias{mix.compnorm}
\title{
Gaussian mixture models for compositional data
}

\description{
Gaussian mixture models for compositional data.
}

\usage{
mix.compnorm(x, g, model, type = "alr", veo = FALSE)
}

\arguments{
\item{x}{
A matrix with the compositional data.
}
\item{g}{
How many clusters to create.
}
\item{model}{
The type of model to be used.
\enumerate{
\item "EII": All groups have the same diagonal covariance matrix, with the same variance for all variables.
\item "VII": Different diagonal covariance matrices, with the same variance for all variables within each group.
\item "EEI": All groups have the same diagonal covariance matrix.
\item "VEI": Different diagonal covariance matrices. If we make all covariance matrices have determinant 1,
(divide the matrix
with the $p$-th root of its determinant) then all covariance matrices will be the same.
\item "EVI": Different diagonal covariance matrices with the same determinant.
\item "VVI": Different diagonal covariance matrices, with nothing in common.
\item "EEE": All covariance matrices are the same.
\item "EEV": Different covariance matrices, but with the same determinant and in addition, if we make them
have determinant 1,
they will have the same trace.
\item "VEV": Different covariance matrices but if we make the matrices have determinant 1, then they will
have the same trace.
\item "VVV": Different covariance matrices with nothing in common.
\item "EVE": Different covariance matrices, but with the same determinant. In addition, calculate the
eigenvectors for each covariance matrix and you will see the extra similarities.
\item "VVE": Different covariance matrices, but they have something in common with their directions.
Calculate the eigenvectors of each covariance matrix and you will see the similarities.
\item "VEE": Different covariance matrices, but if we make the matrices have determinant 1, then they will
have the same trace.
In addition, calculate the eigenvectors for each covariance matrix and you will see the extra similarities.
\item "EVV": Different covariance matrices, but with the same determinant.
}
}
\item{type}{
The type of trasformation to be used, either the additive log-ratio ("alr"), the isometric log-ratio ("ilr")
or the pivot coordinate ("pivot") transformation.
}
\item{veo}{
Stands for "Variables exceed observations". If TRUE then if the number variablesin the model exceeds
the number
of observations, but the model is still fitted.
}
}

\details{
A log-ratio transformation is applied and then a Gaussian mixture model is constructed.
}

\value{
A list including:
\item{mu}{
A matrix where each row corresponds to the mean vector of each cluster.
}
\item{su}{
An array containing the covariance matrix of each cluster.
}
\item{prob}{
The estimated mixing probabilities.
}
\item{est}{
The estimated cluster membership values.
}
}

\references{
Ryan P. Browne, Aisha ElSherbiny and Paul D. McNicholas (2015). R package mixture: Mixture Models
for Clustering and Classification.

Aitchison J. (1986). The statistical analysis of compositional data. Chapman & Hall.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{bic.mixcompnorm}, \link{rmixcomp}, \link{mix.compnorm.contour}, \link{alfa.mix.norm},
\link{alfa.knn},
\link{alfa.rda}, \link{comp.nb}
}
}

\examples{
\dontrun{
x <- as.matrix(iris[, 1:4])
x <- x/ rowSums(x)
mod1 <- mix.compnorm(x, 3, model = "EII" )
mod2 <- mix.compnorm(x, 4, model = "VII")
}
}
