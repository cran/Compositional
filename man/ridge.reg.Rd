\name{Ridge regression}
\alias{ridge.reg}
\title{
  Ridge regression
}

\description{
  Ridge regression.
}

\usage{
ridge.reg(y, x, lambda, B = 1, xnew = NULL)
}

\arguments{
  \item{y}{
    A real valued vector. If it contains percentages, the logit transformation is applied.
  }
  \item{x}{
    A matrix with the predictor variable(s), they have to be continuous.
  }
  \item{lambda}{
    The value of the regularisation parameter \eqn{\lambda}.
  }
  \item{B}{
    If B = 1 (default value) no bootstrpa is performed. Otherwise bootstrap standard errors are returned.
  }
  \item{xnew}{
    If you have new data whose response value you want to predict put it here, otherwise leave it as is.
  }
}

\details{
This is used in the function \code{\link{alfa.ridge}}. There is also a built-in function available from the MASS library, called \code{\link{lm.ridge}}.
}

\value{
  A list including:
    \item{beta}{
      The beta coefficients.
    }
  \item{seb}{
    The standard eror of the coefficiens. If B > 1 the bootstrap standard errors will be returned.
  }
  \item{est}{
    The fitted or the predicted values (if xnew is not NULL).
  }
}

\references{
Hoerl A.E. and R.W. Kennard (1970). Ridge regression: Biased estimation for nonorthogonal problems. Technometrics, 12(1): 55-67.

Brown P. J. (1994). Measurement, Regression and Calibration. Oxford Science Publications.
}

\author{
  Michail Tsagris

  R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <gioathineou@gmail.com>
}

%\note{
  %%  ~~further notes~~
    %}

\seealso{
  \code{\link{ridge.tune}, \link{alfa.ridge}, \link{ridge.plot}
}
}

\examples{
y <- as.vector(iris[, 1])
x <- as.matrix(iris[, 2:4])
mod1 <- ridge.reg(y, x, lambda = 0.1)
mod2 <- ridge.reg(y, x, lambda = 0)
}

\keyword{ ridge regression }
\keyword{ regularisation }
