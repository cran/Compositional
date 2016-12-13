\name{Multivariate analysis of variance (James test)}
\alias{maovjames}
\title{
Multivariate analysis of variance (James test)
}

\description{
Multivariate analysis of variance without assuming equality of the covariance matrices.
}

\usage{
maovjames(x, ina, a = 0.05)
}

\arguments{
  \item{x}{
  A matrix containing Euclidean data.
}
  \item{ina}{
A numerical or factor variable indicating the groups of the data.
}
\item{a}{
The significance level, set to 0.005 by default.
}
}

\details{
Multivariate analysis of variance without assuming equality of the covariance matrices.
}

\value{
A vector with the next 4 elements:
\item{test}{
The test statistic.
}
\item{correction}{
The value of the correction factor.
}
\item{corr.critical}{
The corrected critical value of the chi-square distribution.
}
\item{p-value}{
The p-value of the corrected test statistic.
}
}

\references{
G.S.James (1954). Tests of Linear Hypothese in Univariate and Multivariate Analysis when the Ratios of the Population Variances are Unknown. Biometrika, 41(1/2): 19-43.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{maov}, \link{hotel2T2}, \link{james}, \link{comp.test}
}
}

\examples{
maov( as.matrix(iris[,1:4]), iris[,5] )
maovjames( as.matrix(iris[,1:4]), iris[,5] )
}

\keyword{ Multivariate analysis of variance }
\keyword{ No equality of the covariance matrices }
