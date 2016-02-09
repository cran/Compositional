\name{diri.nr}
\alias{diri.nr}
\title{
Fiiting a Dirichlet distribution via Newton-Rapshon
}

\description{
Fiiting a Dirichlet distribution via Newton-Rapshon.
}

\usage{
diri.nr(x)
}

\arguments{
  \item{x}{
  A matrix containing the compositional data.
}
}

\details{
Maximum likelihood estimation of the parameters of a Dirichlet distribution is performed via Newton-Raphson.
Initial values suggested by Minka (2003) are used. The esimtation is much faster than "diri.est" when the sample size and
or the dimensions increase.
}

\value{
A list including:
\item{loglik}{
The value of the log-likelihood.
}
\item{param}{
The estimated parameters.
}
}

\references{
Thomas P. Minka (2003). Estimating a Dirichlet distribution. http://research.microsoft.com/en-us/um/people/minka/papers/dirichlet/minka-dirichlet.pdf
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{diri.est}, \link{diri.contour} \link{rdiri}
}
}

\examples{
x = rdiri( 100, c(5, 7, 1, 3, 10, 2, 4) )
diri.nr(x)
diri.est(x)
}

\keyword{ Dirichlet distribution }
\keyword{ maximum likelihood estimation }
\keyword{ Newton-Raphson }
