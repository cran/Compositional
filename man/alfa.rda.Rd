\name{Regularised and flexible discriminant analysis for compositional data using the alpha-transformation}
\alias{alfa.rda}
\alias{alfa.fda}
\title{
Regularised and flexible discriminant analysis for compositional data using the \eqn{\alpha}-transformation
}

\description{
Regularised and flexible discriminant analysis for compositional data using the \eqn{\alpha}-transformation.
}

\usage{
alfa.rda(xnew, x, ina, a, gam = 1, del = 0)
alfa.fda(xnew, x, ina, a)
}

\arguments{
\item{xnew}{
A matrix with the new compositional data whose group is to be predicted. Zeros are allowed, but you must be carefull to choose strictly positive vcalues of \eqn{\alpha}.
}
\item{x}{
A matrix with the available compositional data. Zeros are allowed, but you must be carefull to choose strictly positive vcalues of \eqn{\alpha}.
}
\item{ina}{
A group indicator variable for the avaiable data.
}
\item{a}{
The value of \eqn{\alpha} for the \eqn{\alpha}-transformation.
}
\item{gam}{
This is a number between 0 and 1. It is the weight of the pooled covariance and the diagonal matrix.
}
\item{del}{
This is a number between 0 and 1. It is the weight of the LDA and QDA.
}
}

\details{
For the alfa.rda, the covariance matrix of each group is calcualted and then the pooled covariance matrix. The spherical covariance matrix consists of the average of the pooled variances in its diagonal and zeros in the off-diagonal elements. gam is the weight of the pooled covariance matrix and 1-gam is the weight of the spherical covariance matrix, Sa = gam * Sp + (1-gam) * sp. Then it is a compromise between LDA and QDA. del is the weight of Sa and 1-del the weight of each group covariance group. This function is a wrapper for \link{alfa.rda}.

For the alfa.fda a flexible discriminant analysis is performed. See the R package \bold{fda} for more details.
}

\value{
For the alfa.rda a list including:
\item{prob}{
The estimated probabilities of the new data of belonging to each group.
}
\item{scores}{
The estimated socres of the new data of each group.
}
\item{est}{
The estimated group membership of the new data.
}
For the alfa.fda a list including:
\item{mod}{
A fda object as returned by the command fda of the R package mda.
}
\item{est}{
The estimated group membership of the new data.
}
}

\references{
Friedman Jerome, Trevor Hastie and Robert Tibshirani (2009). The elements of statistical learning, 2nd edition. Springer, Berlin.

Tsagris Michail, Simon Preston and Andrew T.A. Wood (2016). Improved classification for compositional data using the \eqn{\alpha}-transformation. Journal of classification, 33(2): 243-261.
https://arxiv.org/pdf/1106.1451.pdf

Hastie, Tibshirani and Buja (1994). Flexible Disriminant Analysis by Optimal Scoring. Journal of the American Statistical Association, 89(428):1255-1270.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} and Giorgos Athineou <gioathineou@gmail.com>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{rda}, \link{alfa}
}
}

\examples{
x <- as.matrix(iris[, 1:4])
x <- x / rowSums(x)
ina <- iris[, 5]
mod <- alfa.rda(x, x, ina, 0)
table(ina, mod$est)
mod2 <- alfa.fda(x, x, ina, 0)
table(ina, mod2$est)
}

\keyword{ Regularised discriminant analysis }

