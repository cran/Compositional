\name{Tuning of the bandwidth h of the kernel using the maximum likelihood cross validation}
\alias{mkde.tune}
\title{
Tuning of the bandwidth h of the kernel using the maximum likelihood cross validation
}

\description{
Tuning of the bandwidth h of the kernel using the maximum likelihood cross validation.
}

\usage{
mkde.tune( x, low = 0.1, up = 3, s = cov(x) )
}

\arguments{
  \item{x}{
  A matrix with Euclidean (continuous) data.
}
  \item{low}{
  The minimum value to search for the optimal bandwidth value.
}
  \item{up}{
  The maximum value to search for the optimal bandwidth value.
}
  \item{s}{
  A covariance matrix. By default it is equal to the covariance matrix of the data, but can change to a robust covariance matrix, MCD for example.
}
}

\details{
Maximum likelihood cross validation is applied in order to choose the optimal value of the bandwidth parameter. No plot is produced.
}

\value{
A list including:
  \item{hopt}{
  The optimal bandwidth value.
}
  \item{maximum}{
  The value of the pseudo-log-likelihood at that given bandwidth value.
}

}
\references{
Arsalane Chouaib Guidoum (2015). Kernel Estimator and Bandwidth Selection for Density and its Derivatives. The kedd package.
http://cran.r-project.org/web/packages/kedd/vignettes/kedd.pdf

M.P. Wand and M.C. Jones (1995). Kernel smoothing, pages 91-92.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <gioathineou@gmail.com>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{mkde}, \link{comp.kerncontour}
}
}

\examples{
library(MASS)
mkde.tune(as.matrix(iris[, 1:4]), c(0.1, 3) )
}

\keyword{ maximum likelihood cross validation }
\keyword{ bandwidth tuning }
\keyword{ kernel density }
