\name{MLE of distributions defined in the (0, 1) interval}
\alias{beta.est}
\alias{ibeta.est}
\alias{logitnorm.est}
\alias{hsecant01.est}
\title{
MLE of distributions defined in the (0, 1) interval
}

\description{
MLE of distributions defined in the (0, 1) interval.
}

\usage{
beta.est(x, tol = 1e-09)
ibeta.est(x, tol = 1e-09)
logitnorm.est(x)
hsecant01.est(x, tol = 1e-09)
}

\arguments{
\item{x}{
A numerical vector with proportions, i.e. numbers in (0, 1) (zeros and ones are not allowed).
}
\item{tol}{
The tolerance level up to which the maximisation stops.
}
}

\details{
Maximum likelihood estimation of the parameters of the beta  distribution is performed via Newton-Raphson. The distributions and hence the functions
does not accept zeros. "logitnorm.mle" fits the logistic normal, hence no nwewton-Raphson is required and the "hypersecant01.mle" uses the golden ratio
search as is it faster than the Newton-Raphson (less calculations)
}

\value{
A list including:
\item{iters}{
The number of iterations required by the Newton-Raphson.
}
\item{loglik}{
The value of the log-likelihood.
}
\item{param}{
The estimated parameters. In the case of "hypersecant01.est" this is called "theta" as there is only one parameter.
}
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Manos Papadakis <papadakm95@gmail.com>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{diri.nr2},
}
}

\examples{
x <- rbeta(1000, 1, 4)
beta.est(x)
ibeta.est(x)

x <- runif(1000)
hsecant01.est(x)
logitnorm.est(x)
ibeta.est(x)

x <- rbeta(1000, 2, 5)
x[sample(1:1000, 50)] <- 0
ibeta.est(x)




}

\keyword{ Beta distribution }
\keyword{ logistic normal distribution }
\keyword{ inflated beta distribution }
\keyword{ hypersecant distribution for proportions }
\keyword{ maximum likelihood estimation }
\keyword{ Newton-Raphson }
\keyword{ proportions }