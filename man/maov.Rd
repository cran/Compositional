\name{Multivariate analysis of variance assuming equality of the covariance matrices}
\alias{maov}
\title{
Multivariate analysis of variance assuming equality of the covariance matrices
}

\description{
Multivariate analysis of variance assuming equality of the covariance matrices.
}

\usage{
maov(x, ina)
}

\arguments{
\item{x}{
A matrix containing Euclidean data.
}
\item{ina}{
A numerical or factor variable indicating the groups of the data.
}
}

\details{
Multivariate analysis of variance assuming equality of the covariance matrices.
}

\value{
A list including:
\item{note}{
A message stating whether the \eqn{F} or the \eqn{chi^2} approximation has been used.
}
\item{result}{
The test statistic and the p-value.
}
}

\references{
Johnson R.A. and Wichern D.W. (2007, 6th Edition). Applied Multivariate Statistical Analysis, pg. 302--303.

Todorov V. and Filzmoser P. (2010). Robust Statistic for the One-way MANOVA.
Computational Statistics & Data Analysis, 54(1): 37--48.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{maovjames}, \link{hotel2T2}, \link{james}
}
}

\examples{
maov( as.matrix(iris[,1:4]), iris[,5] )
maovjames( as.matrix(iris[,1:4]), iris[,5] )
}

