\name{Ternary diagram with the coefficients of the simplicial-simplicial regression models}
\alias{ternary.coef}

\title{
Ternary diagram with the coefficients of the simplicial-simplicial regression models
}

\description{
Ternary diagram with the coefficients of the simplicial-simplicial regression models.
}

\usage{
ternary.coef(B, dg = FALSE, hg = FALSE, colour = NULL)
}

\arguments{
\item{B}{
A matrix with the coefficients of the \code{\link{tflr}} or the \code{\link{scls}}
functions. See examples for this.
}
\item{dg}{
Do you want diagonal grid lines to appear? If yes, set this TRUE.
}
\item{hg}{
Do you want horizontal grid lines to appear? If yes, set this TRUE.
}
\item{colour}{
If you want the points to appear in different colour put a vector with the colour
numbers or colours.
}
}

\details{
Ternary plot of the coefficients of the \code{\link{tflr}} or
the \code{\link{scls}} functions.
}

\value{
A ternary plot of the coefficients of the \code{\link{tflr}} or the \code{\link{scls}} functions.
}

\references{
Aitchison, J. (1983). Principal component analysis of compositional data. Biometrika 70(1): 57--65.

Aitchison J. (1986). The statistical analysis of compositional data. Chapman & Hall.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{ternary}, \link{ternary.reg}, \link{scls}
}
}

\examples{
y <- as.matrix(iris[, 1:3])
y <- y / rowSums(y)
x <- rdiri(150, runif(5, 1,4) )
mod <- scls(y, x)
ternary.coef(mod$be)
}
