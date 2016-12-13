\name{Tuning the parameters of the regularised discriminant analysis}
\alias{rda.tune}
\title{
Tuning the parameters of the regularised discriminant analysis
}

\description{
Tuning the parameters of the regularised discriminant analysis for Eucldiean data.
}

\usage{
rda.tune(x, ina, M = 10, gam = seq(0, 1, by = 0.1), del = seq(0, 1, by = 0.1),
ncores = 1, mat = NULL)
}

\arguments{
  \item{x}{
  A matrix with the data.
}
  \item{ina}{
  A group indicator variable for the avaiable data.
}
  \item{M}{
  The number of folds in the cross validation.
}
  \item{gam}{
  A grid of values for the \eqn{\gamma} parameter as defined in Tsagris et al. (2016).
}
  \item{del}{
  A grid of values for the \eqn{\delta} parameter as defined in Tsagris et al. (2016).
}
  \item{ncores}{
  The number of cores to use. If more than 1, parallel computing will take place. It is advisable to use it if you have many observations and or many variables, otherwise it will slow down th process.
}
  \item{mat}{
You can specify your own folds by giving a mat, where each column is a fold. Each column contains indices of the observations.  You can also leave it NULL and it will create folds.
}
}

\details{
Cross validation is performed to select the optimal parameters for the regularisded discriminant analysis and also estimate the rate of accuracy using the bias correction method of Tibshirani and Tibshirani (2009).

The covariance matrix of each group is calcualted and then the pooled covariance matrix. The spherical covariance matrix consists of the average of the pooled variances in its diagonal and zeros in the off-diagonal elements. gam is the weight of the pooled covariance matrix and 1-gam is the weight of the spherical covariance matrix, Sa = gam * Sp + (1-gam) * sp. Then it is a compromise between LDA and QDA. del is the weight of Sa and 1-del the weight of each group covariance group. This function is a wrapper for \link{alfa.rda}.
}

\value{
A list including:
If graph is TRUE a plot of the performance versus the number of principal components will appear.
  \item{per}{
  An array with the estimate rate of correct classification for every fold. For each of the M matrices, the row values correspond to gam and the columns to the del parameter.
}
  \item{percent}{
  A matrix with the mean estimated rates of correct classification. The row values correspond to gam and the columns to the del parameter.
}
  \item{se}{
  A matrix with the standard error of the mean estimated rates of correct classification. The row values correspond to gam and the columns to the del parameter.
}
  \item{result}{
  The bias corrected estimated rate of correct classification along with the best gam and del parameters.
}
\item{runtime}{
The time required by the cross-validation procedure.
}
}

\references{
Friedman Jerome, Trevor Hastie and Robert Tibshirani (2009). The elements of statistical learning, 2nd edition. Springer, Berlin

Tsagris Michail, Simon Preston and Andrew TA Wood (2016).  Improved classification for compositional data using the \eqn{\alpha}-transformation. Journal of classification (to appear).
http://arxiv.org/pdf/1106.1451.pdf

Tibshirani and Tibshirani (2009). A bias correction for the minimum error rate in cross-validation. The Annals of Applied Statistics, 3(1):822-829.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{rda}, \link{alfa}
}
}

\examples{
mod <- rda.tune(as.matrix(iris[, 1:4]), iris[, 5], M = 10, gam = seq(0, 1, by = 0.2),
del = seq(0, 1, by = 0.2), ncores = 1, mat = NULL)
mod
}

\keyword{ Regularised discriminant analysis }
\keyword{ Tuning of the hyper-parameters }
