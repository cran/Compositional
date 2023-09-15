\name{Projections based test for distributional equality of two groups}
\alias{dptest}
\title{
Projections based test for distributional equality of two groups
}

\description{
Projections based test for distributional equality of two groups.
}

\usage{
dptest(x1, x2, B = 100)
}

\arguments{
\item{x1}{
A matrix containing compositional data of the first group.
}
\item{x2}{
A matrix containing compositional data of the second group.
}
\item{B}{
The number of random uniform projections to use.
}
}

\details{
The test compares the distributions of two compositional datasets 
using random projections. For more details see 
Cuesta-Albertos, Cuevas and Fraiman (2009). 
}

\value{
A vector including:
\item{pvalues}{
The p-values of the Kolmogorov-Smirnov tests.
}
\item{pvalue}{
The p-value of the test based on the Benjamini and Heller (2008) procedure.
}
}

\references{
Cuesta-Albertos J. A., Cuevas A. and Fraiman, R. (2009). 
On projection-based tests for directional and compositional data. 
Statistics and Computing, 19: 367--380.

Benjamini Y. and Heller R. (2008). Screening for partial conjunction hypotheses. 
Biometrics, 64(4): 1215--1222.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{comp.test}
}
}

\examples{
x1 <- rdiri(50, c(3, 4, 5)) ## Fisher distribution with low concentration
x2 <- rdiri(50, c(3, 4, 5))
dptest(x1, x2)
}
