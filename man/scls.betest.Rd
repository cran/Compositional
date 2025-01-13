\name{Permutation test for the matrix of coefficients in the SCLS model}
\alias{scls.betest}
\title{
Permutation test for the matrix of coefficients in the SCLS model
}

\description{
Permutation test for the matrix of coefficients in the SCLS model.
}

\usage{
scls.betest(y, x, B, R = 999)
}

\arguments{
\item{y}{
A matrix with the compositional data (dependent variable). Zero values are allowed.
}
\item{x}{
A matrix with the compositional predictors. Zero values are allowed.
}
\item{B}{
A specific matrix of coefficients to test. Under the null hypothesis, the matrix of coefficients is equal to this matrix.
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
Tsagris. M. (2025). Constrained least squares simplicial-simplicial regression.
Statistics and Computing, 35(27).
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
  %%  ~~further notes~~
    %}

\seealso{
\code{\link{scls}, \link{scls2}, \link{tflr}, \link{scls.indeptest},
\link{tflr.indeptest} }
}

\examples{
y <- rdiri(100, runif(3, 1, 3) )
x <- rdiri(100, runif(3, 1, 3) )
B <- diag(3)
scls.betest(y, x, B = B, R = 99)
}

