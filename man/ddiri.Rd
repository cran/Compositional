\name{Density values of a Dirichlet distribution}
\alias{ddiri}
\title{
Density values of a Dirichlet distribution
}

\description{
Density values of a Dirichlet distribution.
}

\usage{
ddiri(x, a, logged = TRUE)
}

\arguments{
\item{x}{
A matrix containing compositional data. This can be a vector or a matrix with the data.
}
\item{a}{
A vector of parameters. Its length must be equal to the number of components, or columns of the matrix with the compositional data and all values must be greater than zero.
}
\item{logged}{
A boolean variable specifying whether the logarithm of the density values to be returned. It is set to TRUE by default.
}
}

\details{
The density of the Dirichlet distribution for a vector or a matrix of compositional data is returned.
}

\value{
A vector with the density values.
}

\references{
Ng Kai Wang, Guo-Liang Tian and Man-Lai Tang (2011). Dirichlet and related distributions: Theory, methods and applications. John Wiley & Sons.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{diri.nr}, \link{diri.est}, \link{diri.contour}, \link{rdiri}
}
}

\examples{
x <- rdiri( 100, c(5, 7, 1, 3, 10, 2, 4) )
a <- diri.est(x)
f <- ddiri(x, a$param)
sum(f)
a
}

\keyword{ Dirichlet distribution }
\keyword{ density values }
