\name{Regularised discriminant analysis for Euclidean data}
\alias{rda}
\title{
Regularised discriminant analysis for Euclidean data
}

\description{
Regularised discriminant analysis for Euclidean data.
}

\usage{
rda(xnew, x, ina, gam = 1, del = 0)
}

\arguments{
  \item{xnew}{
  A matrix with the new data whose group is to be predicted. They have to be continuous.
}
  \item{x}{
  A matrix with the available data. They have to be continuous.
}
  \item{ina}{
  A group indicator variable for the avaiable data.
}
  \item{gam}{
 This is a number between 0 and 1. It is the weight of the pooled covariance and the diagonal matrix.
}
  \item{del}{
 This is a number between 0 and 1. It is the weight of the LDA and QDA.
}
}

\details{
The covariance matrix of each group is calcualted and then the pooled covariance matrix. The spherical covariance matrix consists of the average of the pooled variances in its diagonal and zeros in the off-diagonal elements. gam is the weight of the pooled covariance matrix and 1-gam is the weight of the spherical covariance matrix, Sa = gam * Sp + (1-gam) * sp. Then it is a compromise between LDA and QDA. del is the weight of Sa and 1-del the weight of each group covariance group. This function is a wrapper for \link{alfa.rda}.
}

\value{
A list including:
\item{prob}{
The estimated probabilities of the new data of belonging to each group.
}
\item{scores}{
The estimated socres of the new data of each group.
}
\item{est}{
The estimated group membership of the new data.
}
}

\references{
Friedman Jerome, Trevor Hastie and Robert Tibshirani (2009). The elements of statistical learning, 2nd edition. Springer, Berlin

Tsagris Michail, Simon Preston and Andrew TA Wood (2016).  Improved classification for compositional data using the \eqn{\alpha}-transformation. Journal of classification, 33(2): 243-261.
http://arxiv.org/pdf/1106.1451.pdf
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{rda.tune}, \link{alfa}
}
}

\examples{
x <- as.matrix(iris[, 1:4])
ina <- iris[, 5]
mod <- rda(x, x, ina)
table(ina, mod$est)
}

\keyword{ Regularised discriminant analysis }

