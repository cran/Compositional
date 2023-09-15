\name{Tuning of the projection pursuit regression with compositional predictor variables using the alpha-transformation}
\alias{alfapprcomp.tune}
\title{
Tuning of the projection pursuit regression with compositional predictor variables using the \eqn{\alpha}-transformation
}

\description{
Tuning of the projection pursuit regression with compositional predictor variables using the \eqn{\alpha}-transformation.
}

\usage{
alfapprcomp.tune(y, x, nfolds = 10, folds = NULL, seed = NULL,
nterms = 1:10, a = seq(-1, 1, by = 0.1), graph = FALSE)
}

\arguments{
\item{y}{
A numerical vector with the continuous variable.
}
\item{x}{
A matrix with the available compositional data. Zeros are allowed.
}
\item{nfolds}{
The number of folds to use.
}
\item{folds}{
If you have the list with the folds supply it here.
}
\item{seed}{
You can specify your own seed number here or leave it NULL.
}
\item{nterms}{
The number of terms to try in the projection pursuit regression.
}
\item{a}{
A vector with the values of \eqn{\alpha} for the \eqn{\alpha}-transformation.
}
\item{graph}{
If graph is TRUE (default value) a filled contour plot will appear.
}
}

\details{
The function performs tuning of the projection pursuit regression algorithm with
compositional predictor variables using the \eqn{\alpha}-transformation.
}

\value{
A list including:
\item{runtime}{
The run time of the cross-validation procedure.
}
\item{mse}{
The mean squared error of prediction for each number of terms.
}
\item{opt.nterms}{
The number of terms corresponding to the minimum mean squared error of prediction.
}
\item{opt.alpha}{
The value of \eqn{\alpha} corresponding to the minimum mean squared error of prediction.
}
\item{performance}{
The minimum mean squared error of prediction.
}
}

\references{
Friedman, J. H. and Stuetzle, W. (1981). Projection pursuit regression. Journal of the American
Statistical Association, 76, 817-823. doi: 10.2307/2287576.

Tsagris M.T., Preston S. and Wood A.T.A. (2011). A data-based power transformation for compositional data.
In Proceedings of the 4th Compositional Data Analysis Workshop, Girona, Spain.
https://arxiv.org/pdf/1106.1451.pdf
}


\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{alfa.pprcomp}, \link{pprcomp.tune}, \link{compppr.tune}
}
}

\examples{
x <- as.matrix(iris[, 2:4])
x <- x / rowSums(x)
y <- iris[, 1]
mod <- alfapprcomp.tune( y, x, a = c(0, 0.5, 1) )
}


