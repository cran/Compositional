\name{Constrained linear least squares for compositional responses and predictors}
\alias{ols.compcomp}
\title{
Constrained linear least squares for compositional responses and predictors
}

\description{
Constrained linear least squares for compositional responses and predictors.
}

\usage{
ols.compcomp(y, x, rs = 5, tol = 1e-4, xnew = NULL)
}

\arguments{
\item{y}{
A matrix with the compositional data (dependent variable). Zero values are allowed.
}
\item{x}{
A matrix with the compositional predictors. Zero values are allowed.
}
\item{rs}{
The number of times to run the constrained optimisation using different random starting values each time.
}
\item{tol}{
The threshold upon which to stop the iterations of the constrained optimisation.
}
\item{xnew}{
If you have new data use it, otherwise leave it NULL.
}
}

\details{
The function performs least squares regression where the beta coefficients are constained to be positive and sum to 1.
We were inspired by the transformation-free linear regression for compositional responses and predictors of Fiksel,
Zeger and Datta (2020).

}

\value{
A list including:
\item{runtime}{
The time required by the regression.
}
\item{mse}{
The mean squared errors.
}
\item{be}{
The beta coefficients.
}
\item{est}{
The fitted of xnew if xnew is not NULL.
}
}

\references{
Jacob Fiksel, Scott Zeger and Abhirup Datta (2020). A transformation-free linear regression for 
compositional outcomes and predictors. https://arxiv.org/pdf/2004.07881.pdf

%Murteira, Jose MR, and Joaquim JS Ramalho 2016. Regression analysis of multivariate fractional data. Econometric Reviews 35(4): 515-552.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
  %%  ~~further notes~~
    %}

\seealso{
\code{\link{cv.olscompcomp}, \link{tflr}, \link{kl.alfapcr}  }
}

\examples{
library(MASS)
set.seed(1234)
y <- rdiri(214, runif(4, 1, 3))
x <- as.matrix(fgl[, 2:9])
x <- x / rowSums(x)
mod <- ols.compcomp(y, x, rs = 1)
mod
}

