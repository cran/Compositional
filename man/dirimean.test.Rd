\name{Log-likelihood ratio test for a Dirichlet mean vector}
\alias{dirimean.test}
\title{
Log-likelihood ratio test for a Dirichlet mean vector
}

\description{
Log-likelihood ratio test for a Dirichlet mean vector.
}

\usage{
dirimean.test(x, a)
}

\arguments{
\item{x}{
A matrix with the compositional data. No zero values are allowed.
}
\item{a}{
A compositional mean vector. The concentration parameter is estimated at first. If the elements do not sum to 1, it is assumed that the Dirichlet parameters are supplied.
}
}

\details{
Log-likelihood ratio test is performed for the hypothesis the given vector of parameters "a" describes the compositional data well.
}

\value{
If there are no zeros in  the data, a list including:
\item{param}{
A matrix with the estimated parameters under the null and the alternative hypothesis.
}
\item{loglik}{
The log-likelihood under the alternative and the null hypothesis.
}
\item{info}{
The value of the test statistic and its relevant p-value.
}
}

\references{
Ng Kai Wang, Guo-Liang Tian and Man-Lai Tang (2011). Dirichlet and related distributions: Theory, methods and applications. John Wiley \& Sons.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <gioathineou@gmail.com>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{sym.test}, \link{diri.nr}, \link{diri.est}, \link{rdiri}, \link{ddiri}
}
}

\examples{
x <- rdiri( 100, c(1, 2, 3) )
dirimean.test(x, c(1, 2, 3) )
dirimean.test( x, c(1, 2, 3)/6 )
}

\keyword{ Dirichelt mean vector }
\keyword{ Log-likelihood ratio test }
