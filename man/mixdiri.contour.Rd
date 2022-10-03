\name{Contour plot of mixtures of Dirichlet distributions in S^2}
\alias{mixdiri.contour}
\title{
Contour plot of mixtures of Dirichlet distributions in \eqn{S^2}
}

\description{
Contour plot of mixtures of Dirichlet distributions in \eqn{S^2}.
}

\usage{
mixdiri.contour(a, prob, n = 100, x = NULL, cont.line = FALSE)
}

\arguments{
\item{a}{
A matrix where each row contains the parameters of each Dirichlet disctribution.
}
\item{prob}{
A vector with the mixing probabilities.
}
\item{n}{
The number of grid points to consider over which the density is calculated.
}
\item{x}{
This is either NULL (no data) or contains a 3 column matrix with compositional data.
}
\item{cont.line}{
Do you want the contour lines to appear? If yes, set this TRUE.
}
}

\details{
The user can plot only the contour lines of a Dirichlet with a given vector of parameters,
or can also add the relevant data should he/she wish to.
}

\value{
A ternary diagram with the points and the Dirichlet contour lines.
}
\references{
Ng Kai Wang, Guo-Liang Tian and Man-Lai Tang (2011). Dirichlet and related distributions:
Theory, methods and applications. John Wiley & Sons.

Aitchison J. (1986). The statistical analysis of compositional data. Chapman & Hall.
}

\author{
Michail Tsagris and Christos Adam.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} and
Christos Adam \email{pada4m4@gmail.com}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{diri.contour}, \link{gendiri.contour}, \link{compnorm.contour},
\link{comp.kerncontour}, \link{mix.compnorm.contour},
\link{diri.nr}, \link{dda}
}
}

\examples{
a <- matrix( c(12, 30, 45, 32, 50, 16), byrow = TRUE,ncol = 3)
prob <- c(0.5, 0.5)
mixdiri.contour(a, prob)
}
