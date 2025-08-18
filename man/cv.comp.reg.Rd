\name{cross-validation for some compositional regression models}
\alias{cv.comp.reg}
\title{
cross-validation for some compositional regression models
}

\description{
cross-validation for some compositional regression models.
}

\usage{
cv.comp.reg(y, x, type = "comp.reg", nfolds = 10, folds = NULL, seed = NULL)
}

\arguments{
\item{y}{
A matrix with compositional data. Zero values are allowed for some regression models.
}
\item{x}{
The predictor variable(s).
}
\item{type}{
This can be one of the following: "comp.reg", "robust", "kl.compreg", "js.compreg", "diri.reg" or "zadr".
}
\item{nfolds}{
The number of folds to be used. This is taken into consideration only if the folds argument is not supplied.
}
\item{folds}{
If you have the list with the folds supply it here. You can also leave it NULL and it will create folds.
}
\item{seed}{
If seed is TRUE the results will always be the same.
}
}

\details{
A k-fold cross-validation for a compositional regression model is performed.
}

\value{
A list including:
\item{runtime}{
The runtime of the cross-validation procedure.
}
\item{kl}{
The Kullback-Leibler divergences for all runs.
}
\item{js}{
The Jensen-Shannon divergences for all runs.
}
\item{perf}{
The average Kullback-Leibler divergence and average Jensen-Shannon divergence.
}
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{comp.reg}, \link{kl.compreg}, \link{compppr.tune}, \link{aknnreg.tune} }
}

\examples{
y <- as.matrix( iris[, 1:3] )
y <- y / rowSums(y)
x <- iris[, 4]
mod <- cv.comp.reg(y, x)
}

