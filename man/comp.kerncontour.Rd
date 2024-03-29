\name{Contour plot of the kernel density estimate in S^2}
\alias{comp.kerncontour}
\title{
Contour plot of the kernel density estimate in \eqn{S^2}
}

\description{
Contour plot of the kernel density estimate in \eqn{S^2}.
}

\usage{
comp.kerncontour(x, type = "alr", n = 50, cont.line = FALSE)
}

\arguments{
\item{x}{
A matrix with the compositional data. It has to be a 3 column matrix.
}
\item{type}{
This is either "alr" or "ilr", corresponding to the additive and the isometric log-ratio transformation respectively.
}
\item{n}{
The number of grid points to consider, over which the density is calculated.
}
\item{cont.line}{
Do you want the contour lines to appear? If yes, set this TRUE.
}
}

\details{
The alr or the ilr transformation are applied to the compositional data. Then, the optimal bandwidth using maximum
likelihood cross-validation is chosen. The multivariate normal kernel density is calculated for a grid of points.
Those points are the points on the 2-dimensional simplex. Finally the contours are plotted.
}
\value{
A ternary diagram with the points and the kernel contour lines.
}
\references{
M.P. Wand and M.C. Jones (1995). Kernel smoothing, CrC Press.

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
\code{\link{diri.contour}, \link{mix.compnorm.contour}, \link{bivt.contour}, \link{compnorm.contour}
}
}

\examples{
x <- as.matrix(iris[, 1:3])
x <- x / rowSums(x)
comp.kerncontour(x, type = "alr", n = 20)
comp.kerncontour(x, type = "ilr", n = 20)
}

\keyword{ kernel density }
\keyword{ contour plot }
