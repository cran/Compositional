\name{Dirichlet discriminant analysis}
\alias{dda}
\title{
Dirichlet discriminant analysis
}

\description{
Dirichlet discriminant analysis.
}

\usage{
dda(xnew, x, ina)
}
\arguments{
\item{xnew}{
A matrix with the new compositional predictor data whose class you want to
predict. Zeros are allowed.
}
\item{x}{
A matrix with the available compositional predictor data. Zeros are allowed.
}
\item{ina}{
A vector of data. The response variable, which is categorical
(factor is acceptable).
}
}

\details{
The funcitons performs maximum likelihood discriminant analysis using
the Dirichlet distribution.
}

\value{
A vector with the estimated group.
}

\references{
Friedman J., Hastie T. and Tibshirani R. (2017). The elements of statistical
learning. New York: Springer.

Thomas P. Minka (2003). Estimating a Dirichlet distribution.
http://research.microsoft.com/en-us/um/people/minka/papers/dirichlet/minka-dirichlet.pdf

Ng Kai Wang, Guo-Liang Tian and Man-Lai Tang (2011). Dirichlet and related
distributions: Theory, methods and applications. John Wiley & Sons.

Aitchison J. (1986). The statistical analysis of compositional data.
Chapman & Hall.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{cv.dda}, \link{comp.nb}, \link{alfa.rda}, \link{alfa.knn},
\link{comp.knn}, \link{mix.compnorm}, \link{diri.reg}, \link{zadr}
}
}

\examples{
x <- Compositional::rdiri(100, runif(5) )
ina <- rbinom(100, 1, 0.5) + 1
mod <- dda(x, x, ina )
}

