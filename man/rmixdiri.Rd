\name{Simulation of compositional data from mixtures of Dirichlet distributions}
\alias{rmixdiri}
\title{
Simulation of compositional data from mixtures of Dirichlet distributions
}

\description{
Simulation of compositional data from mixtures of Dirichlet distributions.
}

\usage{
rmixdiri(n, a, prob)
}

\arguments{
\item{n}{
The sample size.
}
\item{a}{
A matrix where each row contains the parameters of each Dirichlet component.
}
\item{prob}{
A vector with the mixing probabilities.
}
}

\details{
A sample from a Dirichlet mixture model is generated.
}

\value{
A list including:
\item{id}{
A numeric variable indicating the cluster of simulated vector.
}
\item{x}{
A matrix containing the simulated compositional data.
}
}

\references{
Ye X., Yu Y. K. and Altschul S. F. (2011). On the inference of
Dirichlet mixture priors for protein sequence comparison.
Journal of Computational Biology, 18(8), 941-954.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{rmixcomp}, \link{mixdiri.contour},
}
}

\examples{
a <- matrix( c(12, 30, 45, 32, 50, 16), byrow = TRUE,ncol = 3)
prob <- c(0.5, 0.5)
x <- rmixdiri(100, a, prob)
}
