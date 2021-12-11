\name{Non parametric zero replacement strategies}
\alias{zeroreplace}
\title{
Non parametric zero replacement strategies
}

\description{
Non parametric zero replacement strategies.
}

\usage{
zeroreplace(x, a = 0.65, delta = NULL, type = "multiplicative")
}

\arguments{
\item{x}{
A matrix with the compositional data.
}
\item{a}{
The replacement value (\eqn{\delta}) will be "a" times the minimum value observed in the compositional data.
}
\item{delta}{
Unless you specify the replacement value \eqn{\delta} here.
}
\item{type}{
This can be any of "multiplicative", "additive" or "simple". See the references for more details.
}
}

\details{
The "additive" is the zero replacement strategy suggested in Aitchison (1986, pg. 269).
All of the three strategies can be found in Martin-Fernandez et al. (2003).
}

\value{
A matrix with the zero replaced compositional data.
}

\references{
Martin-Fernandez J. A., Barcelo-Vidal C. & Pawlowsky-Glahn, V. (2003).
Dealing with zeros and missing values in compositional data sets using nonparametric imputation.
Mathematical Geology, 35(3): 253-278.

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
\code{\link{perturbation}, \link{alfa} }
}

\examples{
x <- as.matrix(iris[1:20, 1:4])
x <- x/ rowSums(x)
x[ sample(1:20, 4),  sample(1:4, 1) ] <- 0
x <- x / rowSums(x)
zeroreplace(x)
}

