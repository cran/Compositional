\name{Dirichlet random values simulation}
\alias{rdiri}
\title{
Dirichlet random values simulation
}

\description{
Dirichlet random values simulation.
}

\usage{
rdiri(n, a)
}

\arguments{
  \item{n}{
  The sample size, a numerical value.
}
  \item{a}{
  A numerical vector with the parameter values.
}
}

\details{
The algorithm is straightforward, for each vector, independent gamma values are generated and then divided by their total sum.
}

\value{
A matrix with the simulated data.
}

\references{
Ng Kai Wang, Guo-Liang Tian and Man-Lai Tang (2011). Dirichlet and related distributions: Theory, methods and applications. John Wiley & Sons.

Aitchison J. (1986). The statistical analysis of compositional data. Chapman & Hall.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{diri.est}, \link{diri.nr}, \link{diri.contour}
}
}

\examples{
x <- rdiri( 100, c(5, 7, 1, 3, 10, 2, 4) )
diri.est(x)
}

\keyword{ Dirichlet distribution }
\keyword{ random values simulation }
