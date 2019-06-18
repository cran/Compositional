\name{Tuning of the projection pursuit regression for compositional data}
\alias{compppr.tune}
\title{
Tuning of the projection pursuit regression for compositional data
}

\description{
Tuning of the projection pursuit regression for compositional data
In addition, estimation of the rate of correct classification via K-fold cross-validation.
}

\usage{
compppr.tune(y, x, nfolds = 10, folds = NULL, seed = FALSE, nterms = 1:10, 
type = "alr", yb = NULL, B = 1000 )

}

\arguments{
\item{y}{
A matrix with the available compositional data, but zeros are not allowed.
}
\item{x}{
A matrix with the continuous predictor variables.
}
\item{nfolds}{
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
Either "alr" or "ilr" corresponding to the additive or the isometric log-ratio transformation respectively.
}
\item{yb}{
If you have already transformed the data using a log-ratio transformation put it here.
Othewrise leave it NULL.
}
\item{B}{
The number of bootstrap re-samples to use for the unbiased estimation of the performance of the 
projection pursuit regression. If B = 1, no bootstrap is applied.
}
}

\details{
The function performs tuning of the projection pursuit regression algorithm. 
}

\value{
A list including:
\item{kl}{
The average Kullback-Leibler divergence.
}
\item{bc.perf}{
The bootstrap bias corrected average Kullback-Leibler divergence. If no bootstrap was performed this is equal to
the average Kullback-Leibler divergence.
}
\item{runtime}{
The run time of the cross-validation procedure.
}
}

\references{
Friedman, J. H. and Stuetzle, W. (1981). Projection pursuit regression. Journal of the American
Statistical Association, 76, 817-823. doi: 10.2307/2287576.

Tsamardinos I., Greasidou E. and Borboudakis G. (2018). 
Bootstrapping the out-of-sample predictions for efficient and accurate cross-validation. 
Machine Learning 107(12): 1895-1922.  
https://link.springer.com/article/10.1007/s10994-018-5714-4
}


\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <gioathineou@gmail.com>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{comp.ppr}, \link{comp.reg}, \link{alfa}
}
}

\examples{
x <- as.matrix(iris[, 1:4])
x <- x/ rowSums(x)
ina <- iris[, 5]
mod1 <- compknn.tune(x, ina, a = seq(1, 1, by = 0.1) )
mod2 <- alfaknn.tune(x, ina, a = seq(-1, 1, by = 0.1) )
}

\keyword{ Regularised discriminant analysis }

