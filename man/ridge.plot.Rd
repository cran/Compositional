\name{Ridge regression plot}
\alias{ridge.plot}
\title{
Ridge regression plot
}

\description{
A plot of the regularised regression coefficients is shown.

}
\usage{
ridge.plot(y, x, lambda = seq(0, 5, by = 0.1) )
}

\arguments{
  \item{y}{
  A numeric vector containing the values of the target variable. If the values are proportions or percentages,
  i.e. strictly within 0 and 1 they are mapped into R using the logit transformation. In any case, they must be continuous only.
}
  \item{x}{
  A numeric matrix containing the continuous variables. Rows are samples and columns are features.
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
}

\author{
 Michail Tsagris

R implementation and documentation: Giorgos Athineou <gioathineou@gmail.com> and Michail Tsagris <mtsagris@yahoo.gr>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{ridge.reg}, \link{ridge.tune}, \link{alfa.ridge}, \link{alfaridge.plot}
}
}

\examples{
y <- as.vector(iris[, 1])
x <- as.matrix(iris[, 2:4])
ridge.plot(y, x, lambda = seq(0, 2, by = 0.1) )
}

\keyword{ ridge regression }
\keyword{ plot }
