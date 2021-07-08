\name{Principal component analysis}
\alias{logpca}
\title{
Principal component analysis
}

\description{
Principal component analysis.
}

\usage{
logpca(x, center = TRUE, scale = TRUE, k = NULL, vectors = FALSE)
}

\arguments{
\item{x}{
A matrix with the compositional data. Zero values are not allowed.
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
The logarithm is applied to the compositional data and PCA is performed.
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
\code{\link{alfa.pca}, \link{alfa.pcr}, \link{kl.alfapcr}
}
}

\examples{
x <- as.matrix(iris[, 1:4])
x <- x/ rowSums(x)
a <- logpca(x)
}


