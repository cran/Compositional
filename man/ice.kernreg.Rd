\name{ICE plot for univariate kernel regression}
\alias{ice.kernreg}
\title{
ICE plot for univariate kernel regression
}

\description{
ICE plot for univariate kernel regression.
}

\usage{
ice.kernreg(y, x, h, type = "gauss", k = 1, frac = 0.1)
}

\arguments{
\item{y}{
A numerical vector with the response values.
}
\item{x}{
A numerical matrix with the predictor variables.
}
\item{h}{
The bandwidth value to consider.
}
\item{type}{
The type of kernel to use, "gauss" or "laplace".
}
\item{k}{
Which variable to select?.
}
\item{frac}{
Fraction of observations to use. The default value is 0.1.
}
}

\details{
This function implements the Individual Conditional Expecation plots of Goldstein et al. (2015).
See the references for more details.
}

\value{
A graph with several curves. The horizontal axis contains the selected variable,
whereas the vertical axis contains the centered predicted values. The black curves
are the effects for each observation and the blue line is their average effect.
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
\code{ \link{ice.pprcomp}, \link{kernreg.tune}, \link{alfa.pcr}, \link{lc.reg}
}
}

\examples{
x <- as.matrix( iris[, 2:4] )
y <- iris[, 1]
ice <- ice.kernreg(y, x, h = 0.1, k = 1)
}


