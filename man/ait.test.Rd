\name{Aitchison's test for two mean vectors and/or covariance matrices}
\alias{ait.test}
\title{
Aitchison's test for two mean vectors and/or covariance matrices
}

\description{
Aitchison's test for two mean vectors and/or covariance matrices.
}

\usage{
ait.test(x1, x2, type = 1, alpha = 0.05)
}

\arguments{
\item{x1}{
A matrix containing the compositional data of the first sample.
Zeros are not allowed.
}
\item{x2}{
A matrix containing the compositional data of the second sample.
Zeros are not allowed.
}
\item{type}{
The type of hypothesis test to perform.
Type=1 refers to testing the equality of the mean vectors and the covariance matrices.
Type=2 refers to testing the equality of the covariance matrices.
Type=2 refers to testing the equality of the mean vectors.
}
\item{alpha}{
The significance level, set to 0.05 by default.
}
}

\details{
The test is described in Aitchison (2003). See the references for more information.
}

\value{
A vector with the test statistic, the p-value, the critical value and the
degrees of freedom of the chi-square distribution.
}

\references{
John Aitchison (2003). The Statistical Analysis of Compositional Data, p. 153-157.
Blackburn Press.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{comp.test}
}
}

\examples{
x1 <- as.matrix(iris[1:50, 1:4])
x1 <- x1 / rowSums(x1)
x2 <- as.matrix(iris[51:100, 1:4])
x2 <- x2 / rowSums(x2)
ait.test(x1, x2, type = 1)
ait.test(x1, x2, type = 2)
ait.test(x1, x2, type = 3)
}

