\name{MLE for the multivariate t distribution}
\alias{multivt}
\title{
MLE for the multivariate t distribution
}

\description{
MLE of the parameters of a multivariate t distribution.
}

\usage{
multivt(y, plot = FALSE)
}

\arguments{
\item{y}{
A matrix with continuous data.
}
\item{plot}{
If plot is TRUE the value of the maximum log-likelihood as a function of the degres of freedom is presented.
}
}

\details{
The parameters of a multivariate t distribution are estimated. This is used by the functions \code{\link{comp.den}} and \code{\link{bivt.contour}}.
}

\value{
A list including:
\item{center}{
The location estimate.
}
\item{scatter}{
The scatter matrix estimate.
}
\item{df}{
The estimated degrees of freedom.
}
\item{loglik}{
The log-likelihood value.
}
\item{mesos}{
The classical mean vector.
}
\item{covariance}{
The classical covariance matrix.
}
}

\references{
Nadarajah, S. and Kotz, S. (2008). Estimation methods for the multivariate t distribution.
Acta Applicandae Mathematicae, 102(1):99-118.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} and Giorgos Athineou <gioathineou@gmail.com>.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{bivt.contour}, \link{comp.den}
}
}

\examples{
x <- as.matrix(iris[, 1:4])
multivt(x)
}

\keyword{ multivariate t distribution }
