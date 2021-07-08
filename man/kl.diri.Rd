\name{Kullback-Leibler divergence and Bhattacharyya distance between two Dirichlet distributions}
\alias{kl.diri}
\title{
Kullback-Leibler divergence and Bhattacharyya distance between two Dirichlet distributions
}

\description{
Kullback-Leibler divergence and Bhattacharyya distance between two Dirichlet distributions.
}

\usage{
kl.diri(a, b, type = "KL")
}

\arguments{
\item{a}{
A vector with the parameters of the first Dirichlet distribution.
}
\item{b}{
A vector with the parameters of the second Dirichlet distribution.
}
\item{type}{
A variable indicating whether the Kullback-Leibler divergence ("KL") or the Bhattacharyya distance ("bhatt") is to be computed.
}
}

\details{
Note that the order is important in the Kullback-Leibler divergence, since this is asymmetric, but not in the Bhattacharyya distance, since it is a metric.
}

\value{
The value of the Kullback-Leibler divergence or the Bhattacharyya distance.
}

\references{
Ng Kai Wang, Guo-Liang Tian and Man-Lai Tang (2011). Dirichlet and related distributions: Theory, methods and applications. John Wiley & Sons.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} and Giorgos Athineou <gioathineou@gmail.com>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{diri.est}, \link{diri.nr} }
}

\examples{
library(MASS)
a <- runif(10, 0, 20)
b <- runif(10, 1, 10)
kl.diri(a, b)
kl.diri(b, a)
kl.diri(a, b, type = "bhatt")
kl.diri(b, a, type = "bhatt")
}

\keyword{ Kullback-Leibler divergence }
\keyword{ Bhattacharyya distance }
\keyword{ Dirichlet distribution }
