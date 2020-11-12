\name{Tuning of the alpha generalised correlations between two compositional datasets}
\alias{acor.tune}
\title{
Tuning of the alpha generalised correlations between two compositional datasets
}

\description{
Tuning of the alpha generalised correlations between two compositional datasets.
}

\usage{
acor.tune(y, x, a, type = "dcor")
}

\arguments{
\item{y}{
A matrix with the compositional data.
}
\item{x}{
A matrix with the compositional data.
}
\item{a}{
The range of values of the power transformation to search for the optimal one. If zero values are present it has to be greater than 0.
}
\item{type}{
the type of correlation to compute, the distance correlation ("edist"), the canonical correlation type 1 ("cancor1") or the canonical correlation type 2 ("cancor2"). See details for more information.
}
}

\details{
The \eqn{\alpha}-transformation is applied toe ach composition and then the distance correlation or the canonical correlation is computed. If type = "cancor1" the function returns the value of \eqn{\alpha} that maximizes the product of the eigenvalues. If type = "cancor2" the function returns the value of \eqn{\alpha} that maximizes the the largest eigenvalue.
}

\value{
A list including:
\item{alfa}{ The optimal value of \eqn{\alpha}.
}
\item{acor}{
The maximum value of the acor.
}
\item{runtime}{
The runtime of the optimization
}
}

%\references{
%  Tsagris M.T., Preston S. and Wood A.T.A. (2011). A data-based power transformation for compositional data.
%  In Proceedings of the 4th Compositional Data Analysis Workshop, Girona, Spain.
%  https://arxiv.org/pdf/1106.1451.pdf
%}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

\seealso{
\code{\link{acor}, \link{alfa.profile}, \link{alfa}, \link{alfainv}
}
}

\examples{
y <- rdiri(30, runif(3) )
x <- rdiri(30, runif(4) )
acor(y, x, a = 0.4)
}

