\name{Ridge regression with the alpha-transformation plot}
\alias{alfaridge.plot}
\title{
Ridge regression plot
}

\description{
A plot of the regularised regression coefficients is shown.

}
\usage{
alfaridge.plot(y, x, a, lambda = seq(0, 5, by = 0.1) )
}

\arguments{
\item{y}{
A numeric vector containing the values of the target variable. If the values are proportions or percentages, i.e. strictly within 0 and 1 they are mapped into R using
the logit transformation. In any case, they must be continuous only.
}
\item{x}{
A numeric matrix containing the continuous variables.
}
\item{a}{
The value of the \eqn{\alpha}-transformation. It has to be between -1 and 1. If there are zero values in the data, you must use a strictly positive value.
}
\item{lambda}{
A grid of values of the regularisation parameter \eqn{\lambda}.
}
}

\details{
For every value of \eqn{\lambda} the coefficients are obtained. They are plotted versus the \eqn{\lambda} values.
}

\value{
A plot with the values of the coefficients as a function of \eqn{\lambda}.
}

\references{
Hoerl A.E. and R.W. Kennard (1970). Ridge regression: Biased estimation for nonorthogonal problems. Technometrics, 12(1): 55-67.

Brown P. J. (1994). Measurement, Regression and Calibration. Oxford Science Publications.

Tsagris M.T., Preston S. and Wood A.T.A. (2011). A data-based power transformation for compositional data.
In Proceedings of the 4th Compositional Data Analysis Workshop, Girona, Spain.
https://arxiv.org/pdf/1106.1451.pdf
}

\author{
Michail Tsagris.

R implementation and documentation: Giorgos Athineou <gioathineou@gmail.com> and Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
  %%  ~~further notes~~
    %}

\seealso{
\code{ \link{ridge.plot}, \link{alfa.ridge} }
}

\examples{
library(MASS)
y <- as.vector(fgl[, 1])
x <- as.matrix(fgl[, 2:9])
x <- x / rowSums(x)
alfaridge.plot(y, x, a = 0.5, lambda = seq(0, 5, by = 0.1) )
}

\keyword{ ridge regression }
\keyword{ plot }
