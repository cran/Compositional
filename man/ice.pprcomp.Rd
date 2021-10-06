\name{ICE plot for projection pursuit regression with compositional predictor variables}
\alias{ice.pprcomp}
\title{
ICE plot for projection pursuit regression with compositional predictor variables
}

\description{
ICE plot for projection pursuit regression with compositional predictor variables.
}

\usage{
ice.pprcomp(y, x, nterms = 3, k = 1, type = "alr")
}

\arguments{
\item{y}{
A numerical vector with the continuous variable.
}
\item{x}{
A matrix with the compositional data. No zero values are allowed.
}
\item{nterms}{
The number of terms to include in the final model.
}
\item{k}{
Which variable to select?.
}
\item{type}{
Either "alr" or "log" corresponding to the additive log-ratio transformation or the simple
logarithm applied to the compositional data.
}
}

\details{
This function implements the Individual Conditional Expecation plots of Goldstein et al. (2015). See the references for more details.
}

\value{
A graph with several curves. The horizontal axis contains the selected variable, whereas the vertical axis contains the centered predicted values. The black curves are the effects for each observation and whe blue line is their average effect. The turquoise
line is the smmothed curve produced by the loess function of the average effect and the selected variable.
}

\references{
https://christophm.github.io/interpretable-ml-book/ice.html

Goldstein, A., Kapelner, A., Bleich, J. and Pitkin, E. (2015). Peeking inside the black box: Visualizing statistical
learning with plots of individual conditional expectation. Journal of Computational and Graphical Statistics 24(1): 44-65.

Friedman, J. H. and Stuetzle, W. (1981). Projection pursuit regression. Journal of the American
Statistical Association, 76, 817-823. doi: 10.2307/2287576.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{pprcomp}, \link{pprcomp.tune}, \link{alfa.pcr}, \link{lc.reg}, \link{comp.ppr}
}
}

\examples{
\dontrun{
x <- as.matrix( iris[1:50, 2:4] )
x <- x/ rowSums(x)
y <- iris[1:50, 1]
pprcomp(y, x)
ice.pprcomp(y, x)
}
}


