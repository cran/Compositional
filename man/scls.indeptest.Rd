\name{Permutation linear independence test in the SCLS model}
\alias{scls.indeptest}
\title{
Permutation linear independence test in the SCLS model
}

\description{
Permutation linear independence test in the SCLS model.
}

\usage{
scls.indeptest(y, x, R = 999)
}

\arguments{
\item{y}{
A matrix with the compositional data (dependent variable). Zero values are allowed.
}
\item{x}{
A matrix with the compositional predictors. Zero values are allowed.
}
\item{R}{
The number of permutations to perform.
}
}

\details{
Permutation independence test in the constrained linear least squares for compositional
responses and predictors is performed. The observed test statistic is the MSE computed by \code{\link{scls}}. Then, the rows of X are permuted B times and each time the constrained OLS is performed and the MSE is computed. The p-value is then computed in the usual way.
}

\value{
The p-value for the test of independence between Y and X.
}

\references{
Tsagris. M. (2024). Constrained least squares simplicial-simplicial regression.
https://arxiv.org/pdf/2403.19835.pdf
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
  %%  ~~further notes~~
    %}

\seealso{
\code{\link{scls}, \link{scls2}, \link{tflr}, \link{scls.betest} }
}

\examples{
library(MASS)
set.seed(1234)
y <- rdiri(214, runif(4, 1, 3))
x <- as.matrix(fgl[, 2:9])
x <- x / rowSums(x)
scls.indeptest(y, x, R = 99)
}

