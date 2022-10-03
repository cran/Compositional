\name{Density values of a generalised Dirichlet distribution}
\alias{dgendiri}
\title{
Density values of a generalised Dirichlet distribution
}

\description{
Density values of a generalised Dirichlet distribution.
}

\usage{
dgendiri(x, a, b, logged = TRUE)
}

\arguments{
\item{x}{
A matrix containing compositional data. This can be a vector or a matrix with the data.
}
\item{a}{
A numerical vector with the shape parameter values of the Gamma distribution.
}
\item{b}{
A numerical vector with the scale parameter values of the Gamma distribution.
}
\item{logged}{
A boolean variable specifying whether the logarithm of the density values to
be returned. It is set to TRUE by default.
}
}

\details{
The density of the Dirichlet distribution for a vector or a matrix of compositional
data is returned.
}

\value{
A vector with the density values.
}

\references{
Ng Kai Wang, Guo-Liang Tian and Man-Lai Tang (2011).
Dirichlet and related distributions: Theory, methods and applications.
John Wiley & Sons.

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
\code{\link{ddiri}, \link{rgendiri}, \link{diri.est}, \link{diri.contour}, \link{rdiri}, \link{dda}
}
}

\examples{
a <- c(1, 2, 3)
b <- c(2, 3, 4)
x <- rgendiri(100, a, b)
y <- dgendiri(x, a, b)
}

\keyword{ Dirichlet distribution }
\keyword{ density values }
