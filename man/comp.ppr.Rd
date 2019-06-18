\name{Projection pursuit regression for compositional data}
\alias{comp.ppr}

\title{
Projection pursuit regression for compositional data
}

\description{
Projection pursuit regression for compositional data.
}

\usage{
comp.ppr(y, x, nterms = 3, type = "alr", xnew = NULL, yb = NULL )
}

\arguments{
\item{y}{
A matrix with the compositional data.
}
\item{x}{
A matrix with the continuous predictor variables or a data frame including categorical
predictor variables.
}
\item{nterms}{
The number of terms to include in the final model.
}
\item{type}{
Either "alr" or "ilr" corresponding to the additive or the isometric log-ratio transformation respectively.
}
\item{xnew}{
If you have new data use it, otherwise leave it NULL.
}
\item{yb}{
If you have already transformed the data using a log-ratio transformation put it here.
Othewrise leave it NULL.
}
}

\details{
This is the standard projection pursuit. See the built-in funciton "ppr" for more details.
}

\value{
A list includign:
\item{runtime}{
The runtime of the regression.
}
\item{mod}{
The produced model as returned by the function "ppr".
}
\item{est}{
The fitted values of xnew if xnew is not NULL.
}
}

\references{
Friedman, J. H. and Stuetzle, W. (1981). Projection pursuit regression. Journal of the American
Statistical Association, 76, 817-823. doi: 10.2307/2287576.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr>.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{compknn.tune}, \link{rda}, \link{alfa}
}
}

\examples{
x <- as.matrix( iris[, 1:4] )
x <- x/ rowSums(x)
ina <- iris[, 5]
mod <- comp.knn(x, x, ina, a = 1, k = 5)
table(ina, mod)
mod2 <- alfa.knn(x, x, ina, a = 1, k = 5)
table(ina, mod2)

}

\keyword{ k-NN algorithm }

