\name{Minimized Kullback-Leibler divergence between Dirichlet and logistic normal}
\alias{kl.diri.normal}
\title{
Minimized Kullback-Leibler divergence between Dirichlet and logistic normal
}

\description{
Minimized Kullback-Leibler divergence between Dirichlet and logistic
normal distributions.
}

\usage{
kl.diri.normal(a)
}

\arguments{
\item{a}{
A vector with the parameters of the Dirichlet parameters.
}
}

\details{
The function computes the minimized Kullback-Leibler divergence from the
Dirichlet distribution to the logistic normal distribution.
}

\value{
The minimized Kullback-Leibler divergence from the
Dirichlet distribution to the logistic normal distribution.
}

\references{
Aitchison J. (1986). The statistical analysis of compositional data, p. 127. Chapman & Hall.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{diri.nr}, \link{diri.contour}, \link{rdiri}, \link{ddiri}, \link{dda}, \link{diri.reg}
}
}

\examples{
a <- runif(5, 1, 5)
kl.diri.normal(a)
}

\keyword{ Dirichlet distribution }
\keyword{ maximum likelihood estimation }
