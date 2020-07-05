\name{Density of the Flexible Dirichlet distribution}
\alias{fd.density}
\title{
Density of the Flexible Dirichlet distribution
}

\description{
Density of the Flexible Dirichlet distribution
}

\usage{
fd.density(x, alpha, prob, tau)
}

\arguments{
\item{x}{
A vector or a matrix with compositional data.
}
\item{alpha}{
A vector of the non-negative alpha parameters.
}
\item{prob}{
A vector of the clusters' probabilities. It must sum to one.
}
\item{tau}{
The non-negative scalar tau parameter.
}
}

\details{
For more information see the references.
}

\value{
The density value(s).
}

\references{
Ongaro, A. and Migliorati, S. (2013) A generalization of the Dirichlet distribution. Journal of Multivariate Analysis, 114, 412--426.

Migliorati, S., Ongaro, A. and Monti, G. S. (2016) A structured Dirichlet mixture model for compositional data: inferential and applicative issues. Statistics and Computing, 1--21.
}

\author{
Michail Tsagris ported from the R package FlexDir. \email{mtsagris@uoc.gr}.
}

\examples{
alpha <- c(12, 11, 10)
prob <- c(0.25, 0.25, 0.5)
tau <- 8
x <- rfd(20, alpha, prob, tau)
fd.density(x, alpha, prob, tau)
}

\seealso{
\code{ \link{fd.est}, \link{rfd} }
}

