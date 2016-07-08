\name{Cross validation for the ridge regression with compositional data as predictor using the alpha-transformation}
\alias{alfaridge.tune}
\title{
  Cross validation for the ridge regression with compositional data as predictor using the \eqn{\alpha}-transformation
}

\description{
  Cross validation for the ridge regression is performed using the TT estimate of bias (Tibshirani and Tibshirani, 2009).
  There is an option for the GCV criterion which is automatic. The predictor variables are compositional data and the \eqn{\alpha}-transformation is applied first.
}

\usage{
alfaridge.tune(y, x, M = 10, a = seq(-1, 1, by = 0.1),
lambda = seq(0, 2, by = 0.1), mat = NULL, ncores = 1,
graph = TRUE, col.nu = 15)
}

\arguments{
  \item{y}{
    A numeric vector containing the values of the target variable. If the values are proportions or percentages,
    i.e. strictly within 0 and 1 they are mapped into R using the logit transformation.
  }
  \item{x}{
    A numeric matrix containing the compositional data, i.e. the predictor variables. Zero values are allowed.
  }
  \item{M}{
    The number of folds. Set to 10 by default.
  }
  \item{a}{
    A vector with the a grid of values of \eqn{\alpha} to be used.
  }
  \item{lambda}{
    A vector with the a grid of values of \eqn{\lambda} to be used.
  }
  \item{mat}{
   You can specify your own folds by giving a mat, where each column is a fold. Each column contains indices of the observations.     You can also leave it NULL and it will create folds.
  }
  \item{ncores}{
    The number of cores to use. If it is more than 1 parallel computing is performed. It is advisable to use it if you have many observations and or many variables, otherwise it will slow down th process.
  }
  \item{graph}{
    If graph is TRUE (default value) a filled contour plot will appear.
  }
  \item{col.nu}{
   A number parameter for the filled contour plot, taken into account only if graph is TRUE.
  }
}

\details{
  A k-fold cross validation is performed and the estimated performance is bias corrected as suggested by Tibshirani and Tibshirani (2009).
}

\value{
If graph is TRUE a fileld contour a filled contour will appear.
A list including:
\item{mspe}{
The MSPE where rows correspond to the \eqn{\alpha} values and the columns to the number of principal components.
}
\item{best.par}{
The best pair of \eqn{\alpha} and \eqn{\lambda}.
}
\item{performance}{
The minimum mean squared error of prediction (bias corrected) and the estimated bias.
}
\item{runtime}{
The run time of the cross-validation procedure.
}
}


\references{
Hoerl A.E. and R.W. Kennard (1970). Ridge regression: Biased estimation for nonorthogonal problems. Technometrics, 12(1):55-67.

Brown P. J. (1994). Measurement, Regression and Calibration. Oxford Science Publications.

Tibshirani R.J., and Tibshirani R. (2009). A bias correction for the minimum error rate in cross-validation. The Annals of Applied Statistics 3(2): 822-829.

Tsagris M.T., Preston S. and Wood A.T.A. (2011). A data-based power transformation for compositional data.
In Proceedings of the 4th Compositional Data Analysis Workshop, Girona, Spain. http://arxiv.org/pdf/1106.1451.pdf
}

\author{
Michail Tsagris

R implementation and documentation: Giorgos Athineou <athineou@csd.uoc.gr> and Michail Tsagris <mtsagris@yahoo.gr>
}

%\note{
  %%  ~~further notes~~
    %}

\seealso{
  \code{ \link{alfa.ridge}, \link{ridge.tune} }
}

\examples{
library(MASS)
y <- fgl[, 1]
x <- fgl[, 2:9]
alfaridge.tune( y, x, M = 10, a = seq(0.1, 1, by = 0.1),
lambda = seq(0, 1, by = 0.1) )
}

\keyword{ \eqn{\alpha}-transformation }
\keyword{ Ridge regression }
\keyword{ Parameters tuning }
