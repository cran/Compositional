\name{Ternary diagram of regression models}
\alias{ternary.reg}

\title{
Ternary diagram of regression models
}

\description{
Ternary diagram of regression models.
}

\usage{
ternary.reg(y, est, id, labs)
}

\arguments{
\item{y}{
A matrix with the compositional data.
}
\item{est}{
A matrix with all fitted compositional data for all regression models, one under the other.
}
\item{id}{
A vector indicating the regression model of each fitted compositional data set.
}
\item{labs}{
The names of the regression models to appea in the legend.
}
}

\details{
The points first appear on the ternary plot. Then, the fitted compositional data appear with different lines for each regression model.
}

\value{
The ternary plot and lines for the fitted values of each regression model.
}

\references{
Aitchison J. (1986). The statistical analysis of compositional data. Chapman \& Hall.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{ternary}, \link{diri.contour}, \link{comp.kerncontour}
}
}

\examples{
x <- cbind(1, rnorm(50) )
a <- exp( x \%*\% matrix( rnorm(6,0, 0.4), ncol = 3) )
y <- matrix(NA, 50, 3)
for (i in 1:50) y[i, ] <- rdiri(1, a[i, ])
est <- comp.reg(y, x[, -1], xnew = x[, -1])$est
ternary.reg(y, est, id = rep(1, 50), labs = "ALR regression")
}

\keyword{ ternary plot }
\keyword{ visualisation }
