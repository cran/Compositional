\name{Tuning the value of alpha in the alpha-regression}
\alias{alfareg.tune}
\title{
Tuning the value of \eqn{\alpha} in the \eqn{\alpha}-regression
}

\description{
Tuning the value of \eqn{\alpha} in the \eqn{\alpha}-regression.
}

\usage{
alfareg.tune(y, x, a = seq(0.1, 1, by = 0.1), nfolds = 10,
folds = NULL, nc = 1, seed = FALSE, graph = FALSE)
}

\arguments{
\item{y}{
A matrix with compositional data. zero values are allowed.
}
\item{x}{
A matrix with the continuous predictor variables or a data frame including categorical predictor variables.
}
\item{a}{
The value of the power transformation, it has to be between -1 and 1. If zero values are present it has to be greater than 0. If \eqn{\alpha=0} the isometric log-ratio transformation is applied.
}
\item{nfolds}{
The number of folds to split the data.
}
\item{folds}{
If you have the list with the folds supply it here. You can also leave it NULL and it will create folds.
}
\item{nc}{
The number of cores to use. IF you have a multicore computer it is advisable to use more than 1. It makes the procedure faster. It is advisable to use it if you have many observations and or many variables, otherwise it will slow down th process.
}
\item{seed}{
If seed is TRUE the results will always be the same.
}
\item{graph}{
If graph is TRUE a plot of the performance for each fold along the values of \eqn{\alpha} will appear.
}
}

\details{
The \eqn{\alpha}-transformation is applied to the compositional data and the numerical optimisation is performed for the regression, unless \eqn{\alpha=0}, where the coefficients are available in closed form.
}

\value{
A plot of the estimated Kullback-Leibler divergences (multiplied by 2) along the values of \eqn{\alpha} (if graph is set to TRUE).
A list including:
\item{runtime}{
The runtime required by the cross-validation.
}
\item{kula}{
A matrix with twice the Kullback-Leibler divergence of the observed from the fitted values. Each row corresponds to a fold and each column to a value of \eqn{\alpha}. The average over the columns equal the next argument, "kl".
}
\item{kl}{
A vector with twice the Kullback-Leibler divergence of the observed from the fitted values. Every value corresponds to a value of \eqn{\alpha}.
}
  \item{opt}{
The optimal value of \eqn{\alpha}.
}
  \item{value}{
The minimum value of twice the Kullback-Leibler.
}
}

\references{
Tsagris Michail (2015). Regression analysis with compositional data containing zero values. Chilean Journal of Statistics, 6(2): 47-57.
http://arxiv.org/pdf/1508.01913v1.pdf

Tsagris M.T., Preston S. and Wood A.T.A. (2011). A data-based power transformation for compositional data.
In Proceedings of the 4th Compositional Data Analysis Workshop, Girona, Spain.
http://arxiv.org/pdf/1106.1451.pdf
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <gioathineou@gmail.com>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{alfa.reg}, \link{alfa}
}
}

\examples{
library(MASS)
y <- as.matrix(fgl[1:40, 2:4])
y <- y /rowSums(y)
x <- as.vector(fgl[1:40, 1])
mod <- alfareg.tune(y, x, a = c(0.2, 0.35, 0.05), nfolds = 5)
}

