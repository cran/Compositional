\name{Contour plot of a Kent distribution in S^2}
\alias{kent.contour}
\title{
Contour plot of a Kent distribution in \eqn{S^2}
}

\description{
Contour plot of a Kent distribution in \eqn{S^2}.
}

\usage{
kent.contour(G, param, n = 100, x = NULL)
}

\arguments{
\item{G}{
A 3 x 3 matrix whose first column is the mean direction. The second and third columns are the major and minor axes respectively.
}
\item{param}{
A vector with the concentration \eqn{kappa} and ovalness \eqn{\beta} parameters (the \eqn{\psi} parameter has been absorbed inside the matrix G).
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
A ternary diagram with the points and the Dirichlet contour lines.
}
\references{
Graf, M. (2020). Regression for compositions based on a generalization of the Dirichlet distribution. Statistical Methods & Applications, (to appear).
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
G <- matrix( c(0.05713382, 0.96029716, 0.27306608, -0.98809661,
0.01525976, 0.15307588, 0.1428314, -0.2785615, 0.9497382), ncol = 3 )
param <- c(2361.8401338, 1171.3808172, 0.1435577)
kent.contour(G, param)
}

