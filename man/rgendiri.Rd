\name{Generalised Dirichlet random values simulation}
\alias{rgendiri}
\title{
Generalised Dirichlet random values simulation
}

\description{
Generalised Dirichlet random values simulation.
}

\usage{
rgendiri(n, a, b)
}

\arguments{
\item{n}{
The sample size, a numerical value.
}
\item{a}{
A numerical vector with the shape parameter values of the Gamma distribution.
}
\item{b}{
A numerical vector with the scale parameter values of the Gamma distribution.
}
}

\details{
The algorithm is straightforward, for each vector, independent gamma values are generated and
then divided by their total sum. The difference with \code{\link{rdiri}} is that
here the Gamma distributed variables are not equally scaled.
}

\value{
A matrix with the simulated data.
}

\references{
Ng Kai Wang, Guo-Liang Tian and Man-Lai Tang (2011). Dirichlet and related distributions: Theory, methods and applications. John Wiley & Sons.

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
\code{\link{rdiri}, \link{diri.est}, \link{diri.nr}, \link{diri.contour}
}
}

\examples{
a <- c(1, 2, 3)
b <- c(2, 3, 4)
x <- rgendiri(100, a, b)
}

\keyword{ Generalised Dirichlet distribution }
\keyword{ random values simulation }
