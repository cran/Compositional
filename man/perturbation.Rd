\name{Perturbation operation}
\alias{perturbation}
\title{
Perturbation operation
}

\description{
Perturbation operation.
}

\usage{
perturbation(x, y, oper = "+")
}

\arguments{
\item{x}{
A matrix with the compositional data.
}
\item{y}{
Either a matrix with compositional data or a vector with compositional data. In either case, the data may not be compositional data, as long as they non negative.
}
\item{oper}{
For the summation this must be "*" and for the negation it must be "/". According to Aitchison (1986), multiplication is equal to summation in the log-space, and division is equal to negation.
}
}

\details{
This is the perturbation operation defined by Aitchison (1986).
}

\value{
A matrix with the perturbed compositional data.
}

\references{
Aitchison J. (1986). The statistical analysis of compositional data. Chapman & Hall.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{power} }
}

\examples{
x <- as.matrix(iris[1:15, 1:4])
y <- as.matrix(iris[21:35, 1:4])
perturbation(x, y)
perturbation(x, y[1, ])
}

