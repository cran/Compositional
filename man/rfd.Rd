\name{Simulation of compositional data from the Flexible Dirichlet distribution}
\alias{rfd}
\title{
Simulation of compositional data from the Flexible Dirichlet distribution
}

\description{
Simulation of compositional data from the Flexible Dirichlet distribution.
}

\usage{
rfd(n, alpha, prob, tau)
}

\arguments{
\item{n}{
The sample size.
}
\item{alpha}{A vector of the non-negative \eqn{alpha} parameters.
}
\item{prob}{
A vector of the clusters' probabilities that must sum to one.
}
\item{tau}{The positive scalar \eqn{tau} parameter.
}
}

\details{
For more information see the references.
}

\value{
A matrix with compositional data.
}

\references{
Ongaro A. and Migliorati S. (2013). A generalization of the Dirichlet distribution.
Journal of Multivariate Analysis, 114, 412--426.

Migliorati S., Ongaro A. and Monti G. S. (2017). A structured Dirichlet mixture model for compositional data: inferential and applicative issues. Statistics and Computing, 27, 963--983.
}

\author{
Michail Tsagris ported from the R package FlexDir. \email{mtsagris@uoc.gr}.
}

\examples{
alpha <- c(12, 11, 10)
prob <- c(0.25, 0.25, 0.5)
x <- rfd(100, alpha, prob, 7)
}

\seealso{
\code{ \link{fd.est}, \link{fd.density} }
}
