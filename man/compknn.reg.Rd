\name{k-NN regression with compositional data}
\alias{compknn.reg}
\title{
k-NN regression with compositional data
}

\description{
k-NN regression with compositional data.
}

\usage{
compknn.reg(y, x, xnew, k = 2:10, type = "alr", yb = NULL)
}

\arguments{
\item{y}{
A matrix with compositional data. Zero values are not allowed.
}
\item{x}{
The predictor variable(s), they have to be continuous.
}
\item{xnew}{
This is by default set to NULL. If you have new data whose compositional data values you want to predict, put them here.
}
\item{k}{
The number of nearest neighbours to consider. It can be a single number or a vector.
}
\item{type}{
This is either "alr" or "ilr", corresponding to the additive and the isometric log-ratio transformation respectively.
}
\item{yb}{
If you have already transformed the data using the additive log-ratio transformation, plut it here. Othewrise leave it NULL.
This is intended to be used in the function \code{\link{alfareg.tune}} in order to speed up the process.
}
}

\details{
The k-NN regression is applied with responses being compositional data and the predictors being euclidean data.
}

\value{
A list including:
\item{est}{
The fitted values of xnew. If you provide more than one k nearest neighbours this will contain a matrix of the estimated
values for each value of k.
}
}

\references{
Aitchison J. (1986). The statistical analysis of compositional data. Chapman \& Hall.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr>.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{comp.reg}, \link{spatmed.reg}, \link{js.compreg}, \link{comp.knn}
}
}

\examples{
library(MASS)
y <- as.matrix(iris[, 1:3])
y <- y / rowSums(y)
x <- as.vector(iris[, 4])
mod <- compknn.reg(y, x, xnew = x, k = c(2:4))
}

