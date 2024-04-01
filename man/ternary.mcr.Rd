\name{Ternary diagram with confidence region for the mean}
\alias{ternary.mcr}

\title{
Ternary diagram with confidence region for the mean
}

\description{
Ternary diagram with confidence region for the mean.
}

\usage{
ternary.mcr(x, type = "alr", conf = 0.95, dg = FALSE, hg = FALSE, colour = NULL)
}

\arguments{
\item{x}{
A matrix with the compositional data.
}
\item{dg}{
Do you want diagonal grid lines to appear? If yes, set this TRUE.
}
\item{type}{
The type of log-ratio transformation to aply, the "alr" or the "ilr".
}
\item{conf}{
The confidence level, by default this is set to 0.95.
}
\item{hg}{
Do you want horizontal grid lines to appear? If yes, set this TRUE.
}
\item{colour}{
If you want the points to appear in different colour put a vector with the colour
numbers or colours.
}
}

\details{
Ternary plot of compositional data including the log-ratio mean and its confidence region.
The confidence region is based on the Hotelling \eqn{T^2} test statistic of the log-ratio
transformed data.
}

\value{
A ternary plot of compositional data including the log-ratio mean and its confidence region.
}

\references{
Aitchison, J. (1983). Principal component analysis of compositional data. Biometrika 70(1):57-65.

Aitchison J. (1986). The statistical analysis of compositional data. Chapman & Hall.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{ternary}, \link{ternary.reg}, \link{diri.contour}
}
}

\examples{
x <- as.matrix(iris[, 1:3])
x <- x / rowSums(x)
ternary.mcr(x, type = "alr", dg = TRUE, hg = TRUE)
}
