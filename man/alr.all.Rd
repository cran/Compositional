\name{All pairwise additive log-ratio transformations}
\alias{alr.all}
\title{
All pairwise additive log-ratio transformations
}

\description{
All pairwise additive log-ratio transformations.
}

\usage{
alr.all(x)
}

\arguments{
\item{x}{
A numerical matrix with the compositional data.
}
}

\details{
The additive log-ratio transformation with the first component being the commn divisor is applied. Then all the other pairwise log-ratios are computed and added next to each column. For example, divide by the first component, then divide by the second component and so on. This means that no zeros are allowed.
}

\value{
A matrix with all pairwise alr transformed data.
}

\references{
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
\code{\link{alr}, \link{alfa}
}
}

\examples{
x <- as.matrix(iris[, 2:4])
x <- x / rowSums(x)
y <- alr.all(x)
}

