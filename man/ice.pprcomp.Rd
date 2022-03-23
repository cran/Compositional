\name{ICE plot for projection pursuit regression with compositional predictor variables}
\alias{ice.pprcomp}
\title{
ICE plot for projection pursuit regression with compositional predictor variables
}

\description{
ICE plot for projection pursuit regression with compositional predictor variables.
}

\usage{
ice.pprcomp(model, x, k = 1, frac = 0.1, type = "log")
}

\arguments{
\item{model}{
The ppr model, the outcome of the \code{\link{pprcomp}} function.
}
\item{x}{
A matrix with the compositional data. No zero values are allowed.
}
\item{k}{
Which variable to select?.
}
\item{frac}{
Fraction of observations to use. The default value is 0.1.
}
\item{type}{
Either "alr" or "log" corresponding to the additive log-ratio transformation or the simple
logarithm applied to the compositional data.
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

Friedman, J. H. and Stuetzle, W. (1981). Projection pursuit regression.
Journal of the American Statistical Association, 76, 817-823. doi: 10.2307/2287576.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{pprcomp}, \link{pprcomp.tune}, \link{ice.kernreg}, \link{alfa.pcr}, \link{lc.reg}, \link{comp.ppr}
}
}

\examples{
x <- as.matrix( iris[, 2:4] )
x <- x/ rowSums(x)
y <- iris[, 1]
model <- pprcomp(y, x)
ice <- ice.pprcomp(model, x, k = 1)
}


