\name{Fitting a Dirichlet distribution}
\alias{diri.est}
\title{
Fitting a Dirichlet distribution
}

\description{
Estimation of the parameters of a fitted Dirichlet distribution.
}

\usage{
diri.est(x, type = "mle")
}

\arguments{
\item{x}{
A matrix containing compositional data.
}
\item{type}{
If you want to estimate the parameters use type="mle". If you want to estimate the mean vector along with the precision parameter,
the second parametrisation of the Dirichlet, use type="prec".
}
}

\details{
Maximum likelihood estimation of the parameters of a Dirichlet distribution is performed.
}

\value{
A list including:
\item{loglik}{
The value of the log-likelihood.
}
\item{param}{
The estimated parameters.
}
\item{phi}{
The estimated precision parameter, if type = "prec".
}
\item{a}{
The estimated mean vector, if type = "prec".
}
\item{runtime}{
The run time of the maximisation procedure.
}
}

\references{
Ng Kai Wang, Guo-Liang Tian and Man-Lai Tang (2011). Dirichlet and related distributions: Theory, methods and applications. John Wiley & Sons.

Aitchison J. (1986). The statistical analysis of compositional data. Chapman & Hall.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} and Giorgos Athineou <gioathineou@gmail.com>.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{diri.nr}, \link{diri.contour}, \link{rdiri}, \link{ddiri}, \link{dda}, \link{diri.reg}
}
}

\examples{
x <- rdiri( 100, c(5, 7, 1, 3, 10, 2, 4) )
diri.est(x)
diri.est(x, type = "prec")
}

\keyword{ Dirichlet distribution }
\keyword{ maximum likelihood estimation }
