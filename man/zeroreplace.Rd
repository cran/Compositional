\name{Aithison's simple zero replacement strategy}
\alias{zeroreplace}
\title{
Aithison's simple zero replacement strategy
}

\description{
Aithison's simple zero replacement strategy.
}

\usage{
zeroreplace(x, a = 2/3)
}

\arguments{
\item{x}{
A matrix with the compositional data.
}
\item{a}{
The replacement value will be "a" times the minimum value observed in the compositional data.
}
}

\details{
This is the simple zero replacement strategy suggested in Aitchison (1986, pg. 269). 
}

\value{
A matrix with the zero replaced compositional data.
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
\code{\link{perturbation}, \link{alfa} }
}

\examples{
x <- as.matrix(iris[1:20, 1:4])
x <- x/ rowSums(x)
x[ sample(1:20, 4),  sample(1:4, 1) ] <- 0
zeroreplace(x)
}

