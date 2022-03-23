\name{Cross validation for the kernel regression with Euclidean response data}
\alias{kernreg.tune}
\title{
Cross validation for the kernel regression with Euclidean response data
}

\description{
Cross validation for the kernel regression with Euclidean response data.
}

\usage{
kernreg.tune(y, x, h = seq(0.1, 1, length = 10), type = "gauss",
nfolds = 10, folds = NULL, seed = NULL, graph = FALSE, ncores = 1)
}

\arguments{
\item{y}{
A matrix or a vector with the Euclidean response.
}
\item{x}{
A matrix with the available predictor variables.
}
\item{h}{
A vector with the bandwidth value(s) \eqn{h} to consider.
}
\item{type}{
The type of kernel to use, "gauss" or "laplace".
}
\item{nfolds}{
The number of folds. Set to 10 by default.
}
\item{folds}{
If you have the list with the folds supply it here. You can also leave it NULL and it will create folds.
}
\item{seed}{
You can specify your own seed number here or leave it NULL.
}
\item{graph}{
If graph is TRUE (default value) a plot will appear.
}
\item{ncores}{
The number of cores to use. Default value is 1.
}
}

\details{
A k-fold cross validation for the kernel regression with a euclidean response is performed.
}

\value{
A list including:
\item{mspe}{
The mean squared prediction error (MSPE) for each fold and value of \eqn{h}.
}
\item{h}{
The optimal \eqn{h} that leads to the minimum MSPE.
}
\item{performance}{
The minimum MSPE.
}
\item{runtime}{
The runtime of the cross-validation procedure.
}
}

\references{
Wand M. P. and Jones M. C. (1994). Kernel smoothing. CRC press.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{kern.reg}, \link{aknnreg.tune}, \link{aknn.reg} }
}

\examples{
y <- iris[, 1]
x <- iris[, 2:4]
mod <- kernreg.tune(y, x, h = c(0.1, 0.2, 0.3) )
}

