\name{Choose the number of principal components via reconstruction error}
\alias{choose.pc}
\title{
Choose the number of principal components via reconstruction error
}

\description{
Choose the number of principal components via reconstruction error.
}

\usage{
choose.pc(x, graph = TRUE)
}

\arguments{
\item{x}{
A numerical matrix with more rows than columns.
}
\item{graph}{
Should the plot of the PRESS values appear? Default value is TRUE.
}
}

\details{
The functions allows for selecting the number of eigenvectors via the reconstruction error which is computed for all eigenvectors based on SVD.
}

\value{
A list including:
\item{values}{
The eigenvalues of the covariance matrix.
}
\item{cumprop}{
The cumulative proportion of the eigenvalues of the covariance matrix.
}
\item{per}{
The differences in the cumulative proportion of the eigenvalues of the covariance matrix.
}
\item{press}{
The reconstruction error \eqn{\sqrt{\sum_{ij}{(x_{ij}-\hat{x}_{ij})^2}}} for each number of eigenvectors.
}
\item{runtime}{
The runtime of the algorithm.
}
}

\references{
Jolliffe I.T. (2002). Principal Component Analysis.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{pcr}, \link{alfa.pcr}, \link{alfapcr.tune}
}
}

\examples{
library(MASS)
x <- as.matrix(fgl[, 2:9])
a <- choose.pc(x, graph = FALSE)
}

