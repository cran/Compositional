\name{Simplicial constrained median regression for compositional responses and predictors model}
\alias{scrq}
\title{
Simplicial constrained median regression for compositional responses and predictors
}

\description{
Simplicial constrained median regression for compositional responses and predictors.
}

\usage{
scrq(y, x, xnew = NULL)
}

\arguments{
\item{y}{
A matrix with the compositional data (dependent variable). Zero values are allowed.
}
\item{x}{
A matrix with the compositional predictors. Zero values are allowed.
}
\item{xnew}{
If you have new data use it, otherwise leave it NULL.
}
}

\details{
The function performs median regression where the beta coefficients are constained to be positive and sum to 1.
}

\value{
A list including:
\item{mlad}{
The mean absolute deviation.
}
\item{be}{
The beta coefficients.
}
\item{est}{
The fitted of xnew if xnew is not NULL.
}
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
\code{\link{scls}, \link{tflr}, \link{scls.indeptest}  }
}

\examples{
library(MASS)
set.seed(1234)
y <- rdiri(214, runif(4, 1, 3))
x <- as.matrix(fgl[, 2:9])
x <- x / rowSums(x)
mod <- scrq(y, x)
mod
}

