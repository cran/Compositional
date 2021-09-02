\name{Fitting a Flexible Dirichlet distribution}
\alias{fd.est}
\title{
Fitting a Flexible Dirichlet distribution
}

\description{
Fitting a Flexible Dirichlet distribution.
}

\usage{
fd.est(x, ini.iter = 50, final.iter = 100)
}

\arguments{
\item{x}{
A matrix or a dataframe containing the compositional data.
}
\item{ini.iter}{
Number of iterations for the initial SEM step. Default value is 50.
}
\item{final.iter}{
Number of iterations for the final EM step. Default value is 100.
}
}

\details{
For more information see the references.
}

\value{
A list including:
\item{alpha}{
Estimated values of the parameter vector \eqn{\alpha}.
}
\item{prob}{
Estimated values of the parameter vector p.
}
\item{tau}{
Estimated value of the parameter \eqn{tau}.
}
\item{loglik}{
The Log-likelihood value.
}
}

\references{
Ongaro A. and Migliorati S. (2013). A generalization of the Dirichlet distribution.
Journal of Multivariate Analysis, 114, 412--426.

Migliorati S., Ongaro A. and Monti G. S. (2017). A structured Dirichlet mixture model for compositional data: inferential and applicative issues. Statistics and Computing, 27, 963--983.
}

\author{
Michail Tsagris ported from the R package FlexDir. \email{mtsagris@uoc.gr}.
}

\examples{
\dontrun{
x <- rfd(n = 50, a = c(12, 11, 10), p = c(0.25, 0.25, 0.5), tau = 4 )
ela <- fd.est(x, ini.iter = 10, final.iter = 20)
ela
}
}

\seealso{
\code{ \link{rfd}, \link{fd.density} }
}
