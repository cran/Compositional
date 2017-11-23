\name{Cross validation for the ridge regression}
\alias{ridge.tune}
\title{
Cross validation for the ridge regression
}

\description{
Cross validation for the ridge regression is performed. There is an option for the GCV criterion which is automatic.
}

\usage{
ridge.tune(y, x, M = 10, lambda = seq(0, 2, by = 0.1), mat = NULL,
ncores = 1, graph = FALSE)
}

\arguments{
  \item{y}{
  A numeric vector containing the values of the target variable. If the values are proportions or percentages,
  i.e. strictly within 0 and 1 they are mapped into R using the logit transformation.
}
  \item{x}{
  A numeric matrix containing the variables.
}
  \item{M}{
  The number of folds. Set to 10 by default.
}
  \item{lambda}{
  A vector with the a grid of values of \eqn{\lambda} to be used.
}
  \item{mat}{
You can specify your own folds by giving a mat, where each column is a fold. Each column contains indices of the observations.  You can also leave it NULL and it will create folds.
}
 \item{ncores}{
  The number of cores to use. If it is more than 1 parallel computing is performed.
}
 \item{graph}{
If graph is set to TRUE the performances for each fold as a function of the \eqn{\lambda} values will appear.
}
}

\details{
A k-fold cross validation is performed. This function is used by \code{\link{alfaridge.tune}}.
}

\value{
A list including:
\item{msp}{
The performance of the ridge regression for every fold.
}
\item{mspe}{
The values of the mean prediction error for each value of \eqn{\lambda}.
}
\item{lambda}{
The value of \eqn{\lambda} which corresponds to the minimum MSPE.
}
\item{performance}{
The minimum MSPE.
}
\item{runtime}{
The time required by the cross-validation procedure.
}
}

\references{
Hoerl A.E. and R.W. Kennard (1970). Ridge regression: Biased estimation for nonorthogonal problems. Technometrics, 12(1):55-67.

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
\code{ \link{ridge.reg}, \link{alfaridge.tune} }
}

\examples{
y <- as.vector(iris[, 1])
x <- as.matrix(iris[, 2:4])
ridge.tune( y, x, M = 10, lambda = seq(0, 2, by = 0.1), graph = TRUE )
}

\keyword{ Ridge regression }
