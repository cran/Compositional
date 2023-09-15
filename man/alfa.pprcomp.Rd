\name{Projection pursuit regression with compositional predictor variables using the alpha-transformation}
\alias{alfa.pprcomp}
\title{
Projection pursuit regression with compositional predictor variables using the \eqn{\alpha}-transformation
}

\description{
Projection pursuit regression with compositional predictor variables using the \eqn{\alpha}-transformation.
}

\usage{
alfa.pprcomp(y, x, nterms = 3, a, xnew = NULL)
}

\arguments{
\item{y}{
A numerical vector with the continuous variable.
}
\item{x}{
A matrix with the compositional data. Zero values are allowed.
}
\item{nterms}{
The number of terms to include in the final model.
}
\item{a}{
The value of \eqn{\alpha} for the \eqn{\alpha}-transformation.
}
\item{xnew}{
If you have new data use it, otherwise leave it NULL.
}
}

\details{
This is the standard projection pursuit. See the built-in function "ppr" for
more details. The compositional data are transformed with the \eqn{\alpha}-transformation
}

\value{
A list including:
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
\code{ \link{alfapprcomp.tune}, \link{pprcomp}, \link{comp.ppr}
}
}

\examples{
x <- as.matrix( iris[, 2:4] )
x <- x / rowSums(x)
y <- iris[, 1]
alfa.pprcomp(y, x, a = 0.5)
}


