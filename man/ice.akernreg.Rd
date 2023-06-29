\name{ICE plot for the alpha-kernel regression}
\alias{ice.akernreg}
\title{
ICE plot for the \eqn{\alpha}-kernel regression
}

\description{
ICE plot for the \eqn{\alpha}-kernel regression.
}

\usage{
ice.akernreg(y, x, a, h, type = "gauss", ind = 1, frac = 0.1, qpos = 0.9)
}

\arguments{
\item{y}{
A numerical vector with the response values.
}
\item{x}{
A numerical matrix with the predictor variables.
}
\item{a}{
The value \eqn{\alpha} to consider.
}
\item{h}{
The bandwidth value to consider.
}
\item{type}{
The type of kernel to use, "gauss" or "laplace".
}
\item{ind}{
Which variable to select?.
}
\item{frac}{
Fraction of observations to use. The default value is 0.1.
}
\item{qpos}{
A number between 0.8 and 1. This is used to place the legend of the figure better.
You can play with it. In the worst case scenario the code is open and you
tweak this argument as you prefer.
}
}

\details{
This function implements the Individual Conditional Expecation plots of Goldstein et al. (2015).
See the references for more details.
}

\value{
A graph with several curves, one for each component. The horizontal axis contains the selected variable,
whereas the vertical axis contains the locally smoothed predicted compositional lines.
}

\references{
https://christophm.github.io/interpretable-ml-book/ice.html

Goldstein, A., Kapelner, A., Bleich, J. and Pitkin, E. (2015).
Peeking inside the black box: Visualizing statistical
learning with plots of individual conditional expectation.
Journal of Computational and Graphical Statistics 24(1): 44-65.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{ice.aknnreg}, \link{ice.pprcomp}
}
}

\examples{
y <- as.matrix( iris[, 2:4] )
x <- iris[, 1]
ice <- ice.akernreg(y, x, a = 0.6, h = 0.1, ind = 1)
}


