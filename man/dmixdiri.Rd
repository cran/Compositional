\name{Density values of a mixture of Dirichlet distributions}
\alias{dmixdiri}
\title{
Density values of a mixture of Dirichlet distributions
}

\description{
Density values of a mixture of Dirichlet distributions.
}

\usage{
dmixdiri(x, a, prob, logged = TRUE)
}

\arguments{
\item{x}{
A vector or a matrix with compositional data. Zeros are not allowed.
}
\item{a}{
A matrix where each row contains the parameters of each Dirichlet component.
}
\item{prob}{
A vector with the mixing probabilities.
}
\item{logged}{
A boolean variable specifying whether the logarithm of the density values to
be returned. It is set to TRUE by default.
}
}

\details{
The density of the mixture of Dirichlet distribution for a vector or a matrix
of compositional data is returned.
}


\value{
A vector with the density values.
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
\code{ \link{rmixdiri}, \link{mixdiri.contour}
}
}

\examples{
a <- matrix( c(12, 30, 45, 32, 50, 16), byrow = TRUE,ncol = 3)
prob <- c(0.5, 0.5)
x <- rmixdiri(100, a, prob)$x
f <- dmixdiri(x, a, prob)
}
