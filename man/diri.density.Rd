\name{Density of the Dirichlet distribution}
\alias{diri.density}
\title{
Density of the Dirichlet distribution
}

\description{
Density of the Dirichlet distribution.
}

\usage{
diri.density(x, a, logged = FALSE)
}

\arguments{
\item{x}{
A vector or a matrix with compositional data.
}
\item{a}{
A vector of the non-negative alpha parameterss.
}
\item{logged}{
Do you want the logarithm of the density values? TRUE or FALSE.
}
}

\value{
The density value(s).
}

%\references{
%}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

\examples{
a <- runif(3, 1, 5)
x <- rdiri(100, a)
a <- diri.nr(x)$param
x <- diri.density(x, a)
}

\seealso{
\code{ \link{fd.density}, \link{rdiri}, \link{diri.nr} }
}

