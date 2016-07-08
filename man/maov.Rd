\name{Multivariate analysis of variance}
\alias{maov}
\title{
Multivariate analysis of variance
}

\description{
Multivariate analysis of variance with assuming equality of the covariance matrices.
}

\usage{
maov(x, ina)
}

\arguments{
  \item{x}{
  A matrix containing the Euclidean data.
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
A message stating whether the F or the chi-sqausre approximation has been used.
}
\item{result}{
The test statistic and the p-value.
}
}

\references{
Johnson and Wichern (2007, 6th Edition). Applied Multivariate Statistical Analysis p. 302-303.

Todorov V. and Filzmoser P. (2010). Robust Statistic for the One-way MANOVA.
Computational Statistics \& Data Analysis 54(1):37-48.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{maovjames}, \link{hotel2T2}, \link{james}, \link{comp.test}
}
}

\examples{
maov( iris[,1:4], iris[,5] )
maovjames( iris[,1:4], iris[,5] )
}

\keyword{ Multivariate analysis of variance }
\keyword{ Equality of the covariance matrices }
