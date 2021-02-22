\name{Principal component analysis using the alpha-transformation}
\alias{alfa.pca}
\title{
Principal component analysis using the \eqn{\alpha}-transformation
}

\description{
Principal component analysis using the \eqn{\alpha}-transformation.
}

\usage{
alfa.pca(x, a, center = TRUE, scale = TRUE, k = NULL, vectors = FALSE)
}

\arguments{
\item{x}{
A matrix with the compositional data. Zero values are allowed. In that case "a" should be positive.
}
\item{a}{
The value of \eqn{\alpha} to use in the \eqn{\alpha}-transformation.
}
\item{center}{
Do you want your data centered? TRUE or FALSE.
}
\item{scale}{
Do you want each of your variables scaled, i.e. to have unit variance? TRUE or FALSE.
}
\item{k}{
If you want a specific number of eigenvalues and eigenvectors set it here, otherwise all
eigenvalues (and eigenvectors if requested) will be returned.
}
\item{vectors}{
Do you want the eigenvectors be returned? By dafault this is FALSE.
}
}

\details{
The \eqn{\alpha}-transformation is applied to the compositional data and then PCA is performed. Note however, that the right multiplication by the Helmert sub-matrix is not applied in order to be in accordance with Aitchison (1983). When \eqn{\alpha=0}, this results to the PCA proposed by Aitchison (1983).
}

\value{
A list including:
\item{values}{
The eigenvalues.
}
\item{vectors}{
The eigenvectors.
}
}

\references{
Aitchison J. (1986). The statistical analysis of compositional data. Chapman \& Hall.

Aitchison, J. (1983). Principal component analysis of compositional data. Biometrika, 70(1), 57-65.

Tsagris M.T., Preston S. and Wood A.T.A. (2011). A data-based power transformation for compositional data.
In Proceedings of the 4th Compositional Data Analysis Workshop, Girona, Spain.
http://arxiv.org/pdf/1106.1451.pdf
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{logpca}, \link{alfa.pcr}, \link{kl.alfapcr}
}
}

\examples{
x <- as.matrix(iris[, 1:4])
x <- x/ rowSums(x)
a <- alfa.pca(x, 0.5)
}


