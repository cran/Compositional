\name{sym.test}
\alias{sym.test}
\title{
Log-likelihood ratio test for a symmetric Dirichlet distribution
}

\description{
Log-likelihood ratio test for a symmetric Dirichlet distribution
}

\usage{
sym.test(x)
}

\arguments{
  \item{x}{
  A matrix with the compositional data. No zero values are allowed.
  }
}

\details{
Log-likelihood ratio test is performed for the hypothesis that all Dirichelt parameters are equal.
}

\value{
A list including:
  \item{est.par}{
The estimated parameters under the alternative hypothesis.
}
  \item{one.par}{
The value of the estimated parameter under the null hypothesis.
}
  \item{res}{
The loglikelihood under the alternative and the null hypothesis, the value of the test statistic, its relevant p-value and the
associated degrees of freedom, which are actually the dimensionality of the simplex, \eqn{D-1}, where \eqn{D} is the number of
components.
}
}

\references{
Ng Kai Wang, Guo-Liang Tian and Man-Lai Tang (2011). Dirichlet and related distributions: Theory, methods and applications. John Wiley \& Sons.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{diri.reg}, \link{esov.compreg}, \link{kl.compreg}, \link{ols.compreg}, \link{comp.reg}
}
}

\examples{
x <- rdiri( 100, c(5, 7, 1, 3, 10, 2, 4) )
sym.test(x)
x <- rdiri( 100, c(5, 5, 5, 5, 5) )
sym.test(x)
}

\keyword{ symmetric Dirichlet distribution }
\keyword{ log-likeliihood ratio test }
