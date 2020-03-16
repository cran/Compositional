\name{Multivariate regression with compositional data}
\alias{comp.reg}
\title{
Multivariate regression with compositional data
}

\description{
Multivariate regression with compositional data.
}

\usage{
comp.reg(y, x, type = "classical", xnew = NULL, yb = NULL)
}

\arguments{
\item{y}{
A matrix with compsitional data. Zero values are not allowed.
}
\item{x}{
The predictor variable(s), they have to be continuous.
}
\item{type}{
The type of regression to be used, "classical" for standard multivariate regression, or "spatial" for the robust spatial median regression. Alternatively you can type "lmfit" for the fast classical multivariate regression that does not return standard errors whatsoever.
}
\item{xnew}{
This is by default set to NULL. If you have new data whose compositional data values you want to predict, put them here.
}
\item{yb}{
If you have already transformed the data using the additive log-ratio transformation, plut it here. Othewrise leave it NULL.
This is intended to be used in the function \code{\link{alfareg.tune}} in order to speed up the process.
}
}

\details{
The additive log-ratio transformation is applied and then the chosen multivariate regression is implemented. The alr is easier to explain than the ilr and that is why the latter is avoided here.
}

\value{
A list including:
\item{runtime}{
The time required by the regression.
}
\item{be}{
The beta coefficients.
}
\item{seb}{
The standard error of the beta coefficients.
}
\item{est}{
The fitted values of xnew if xnew is not NULL.
}
}

\references{
Mardia K.V., Kent J.T., and Bibby J.M. (1979). Multivariate analysis. Academic press.

Aitchison J. (1986). The statistical analysis of compositional data. Chapman \& Hall.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} and Giorgos Athineou <gioathineou@gmail.com>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{multivreg}, \link{spatmed.reg}, \link{js.compreg}, \link{diri.reg}
}
}

\examples{
library(MASS)
y <- as.matrix(iris[, 1:3])
y <- y / rowSums(y)
x <- as.vector(iris[, 4])
mod1 <- comp.reg(y, x)
mod2 <- comp.reg(y, x, type = "spatial")
}

\keyword{ multivariate rgression }
\keyword{ spatial median regression }
