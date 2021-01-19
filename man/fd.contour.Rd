\name{Contour plot of a Flexible Dirichlet distribution in S^2}
\alias{fd.contour}
\title{
Contour plot of a Flexible Dirichlet distribution in \eqn{S^2}
}

\description{
Contour plot of a Flexible Dirichlet distribution in \eqn{S^2}.
}

\usage{
fd.contour(alpha, prob, tau, n = 100, x = NULL)
}

\arguments{
\item{alpha}{
A vector of the non-negative \eqn{alpha} parameters.
}
\item{prob}{
A vector of the clusters' probabilities. It must sum to one.
}
\item{tau}{
The non-negative scalar \eqn{tau} parameter.
}
\item{n}{
The number of grid points to consider over which the density is calculated.
}
\item{x}{
This is either NULL (no data) or contains a 3 column matrix with compositional data.
}
}

\details{
The user can plot only the contour lines of a Dirichlet with a given vector of parameters, or can also add the relevant data should they wish to.
}

\value{
A ternary diagram with the points and the Flexible Dirichlet contour lines.
}

\references{
Ongaro A. and Migliorati S. (2013). A generalization of the Dirichlet distribution.
Journal of Multivariate Analysis, 114, 412--426.

Migliorati S., Ongaro A. and Monti G. S. (2017). A structured Dirichlet mixture model for compositional data: inferential and applicative issues. Statistics and Computing, 27, 963--983.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{norm.contour}, \link{bivt.contour}, \link{comp.kerncontour}, \link{mixnorm.contour}
}
}

\examples{
\dontrun{
x <- rfd(n = 50, a = c(12, 11, 10), p = c(0.25, 0.25, 0.5), tau = 4 )
mod <- fd.est(x, ini.iter = 10, final.iter = 20)
fd.contour(mod$alpha, mod$prob, mod$tau, x = x)
}
}

\keyword{ Dirichlet distribution }
\keyword{ contour plot }
