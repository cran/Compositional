\name{Fitting a Dirichlet distribution via Newton-Rapshon}
\alias{diri.nr}
\title{
Fitting a Dirichlet distribution via Newton-Rapshon
}

\description{
Fitting a Dirichlet distribution via Newton-Rapshon.
}

\usage{
diri.nr(x, type = 1, tol = 1e-07)
}

\arguments{
\item{x}{
A matrix containing compositional data. Zeros are not allowed.
}
\item{type}{
Type can either be 1, so that the Newton-Rapshon is used for the maximisation of the log-likelihood, as Minka (2012) suggested or it 
can be 1. In the latter case the Newton-Raphson algorithm is implemented involving matrix inversions. In addition an even faster 
implementation has been implemented (in C++) in the package \bold{Rfast} and is used here.
}
\item{tol}{
The tolerance level indicating no further increase in the log-likelihood.
}
}

\details{
Maximum likelihood estimation of the parameters of a Dirichlet distribution is performed via Newton-Raphson.
Initial values suggested by Minka (2003) are used. The estimation is super faster than "diri.est" and the
difference becomes really apparent when the sample size and or the dimensions increase. In fact this will work with millions of observations.
So in general, I trust this one more than "diri.est".

The only problem I have seen with this method is that if the data are concentrated around a point,
say the center of the simplex, it will be hard for this and the previous methods to give estimates of the parameters.
In this extremely difficult scenario I would suggest the use of the previous function with the precision parametrization
"diri.est(x, type = "prec")". It will be extremely fast and accurate.
}

\value{
A list including:
\item{iter}{
The number of iterations required. If the argument "type" is set to 2 this is not returned.
}
\item{loglik}{
The value of the log-likelihood.
}
\item{param}{
The estimated parameters.
}
\item{runtime}{
The run time of the procedure.
}
}

\references{
Thomas P. Minka (2003). Estimating a Dirichlet distribution. 
http://research.microsoft.com/en-us/um/people/minka/papers/dirichlet/minka-dirichlet.pdf
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <gioathineou@gmail.com>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{diri.est}, \link{diri.contour} \link{rdiri}, \link{ddiri}
}
}

\examples{
x <- rdiri( 100, c(5, 7, 5, 8, 10, 6, 4) )
diri.nr(x)
diri.nr(x, type = 2)
diri.est(x)
}

\keyword{ Dirichlet distribution }
\keyword{ maximum likelihood estimation }
\keyword{ Newton-Raphson }
