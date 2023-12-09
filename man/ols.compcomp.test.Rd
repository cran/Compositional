\name{Permutation independence test in the constrained linear least squares for compositional responses and predictors}
\alias{ols.compcomp.test}
\title{
Permutation independence test in the constrained linear least squares for compositional responses and predictors
}

\description{
Permutation independence test in the constrained linear least squares for compositional responses and predictors.
}

\usage{
ols.compcomp.test(y, x, B = 999)
}

\arguments{
\item{y}{
A matrix with the compositional data (dependent variable). Zero values are allowed.
}
\item{x}{
A matrix with the compositional predictors. Zero values are allowed.
}
\item{B}{
The number of permutations to perform.
}
}

\details{
Permutation independence test in the constrained linear least squares for compositional
responses and predictors is performed. The observed test statistic is the MSE computed by \code{\link{ols.compcomp}}. Then, the rows of X are permuted B times and each time the constrained OLS is performed and the
MSE is computed. The p-value is then computed in the usual way.
}

\value{
The p-value for the test of independence between Y and X.
}

%\references{
%Jacob Fiksel, Scott Zeger and Abhirup Datta (2020). A transformation-free linear regression for
%compositional outcomes and predictors. https://arxiv.org/pdf/2004.07881.pdf

%Murteira, Jose MR, and Joaquim JS Ramalho 2016. Regression analysis of multivariate fractional data. Econometric Reviews 35(4): 515-552.
%}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
  %%  ~~further notes~~
    %}

\seealso{
\code{\link{ols.compcomp}, \link{tflr} }
}

\examples{
library(MASS)
set.seed(1234)
y <- rdiri(214, runif(4, 1, 3))
x <- as.matrix(fgl[, 2:9])
x <- x / rowSums(x)
ols.compcomp.test(y, x, B = 99)
}

