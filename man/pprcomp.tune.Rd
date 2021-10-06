\name{Tuning of the projection pursuit regression with compositional predictor variables}
\alias{pprcomp.tune}
\title{
Tuning of the projection pursuit regression with compositional predictor variables
}

\description{
Tuning of the projection pursuit regression with compositional predictor variables.
}

\usage{
pprcomp.tune(y, x, nfolds = 10, folds = NULL, seed = FALSE, nterms = 1:10, type = "alr")
}

\arguments{
\item{y}{
A numerical vector with the continuous variable.
}
\item{x}{
A matrix with the available compositional data, but zeros are not allowed.
}\item{nfolds}{
The number of folds to use.
}
\item{folds}{
If you have the list with the folds supply it here.
}
\item{seed}{
If seed is TRUE the results will always be the same.
}
\item{nterms}{
The number of terms to try in the projection pursuit regression.
}
\item{type}{
Either "alr" or "log" corresponding to the additive log-ratio transformation or the logarithm applied to the compositional predictor variables.
}
}

\details{
The function performs tuning of the projection pursuit regression algorithm with compositional predictor variables.
}

\value{
A list including:
\item{runtime}{
The run time of the cross-validation procedure.
}
\item{mse}{
The mean squared error of prediction for each fold and number of terms.
}
}

\references{
Friedman, J. H. and Stuetzle, W. (1981). Projection pursuit regression. Journal of the American
Statistical Association, 76, 817-823. doi: 10.2307/2287576.
}


\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{pprcomp}, \link{ice.pprcomp}, \link{alfapcr.tune}, \link{compppr.tune}
}
}

\examples{
x <- as.matrix(iris[, 2:4])
x <- x/ rowSums(x)
y <- iris[, 1]
mod <- pprcomp.tune(y, x)
}


