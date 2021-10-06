\name{Unconstrained log-contrast regression with compositional predictor variables}
\alias{ulc.reg}
\title{
Unconstrained log-contrast regression with compositional predictor variables
}

\description{
Unconstrained log-contrast regression with compositional predictor variables.
}

\usage{
ulc.reg(y, x, z = NULL, xnew = NULL, znew = NULL)
}

\arguments{
\item{y}{
A numerical vector containing the response variable values. This must be a continuous variable.
}
\item{x}{
A matrix with the predictor variables, the compositional data. No zero values are allowed.
}
\item{z}{
A matrix, data.frame, factor or a vector with some other covariate(s).
}
\item{xnew}{
A matrix containing the new compositional data whose response is to be predicted.
If you have no new data, leave this NULL as is by default.
}
\item{znew}{
A matrix, data.frame, factor or a vector with the values of some other covariate(s).
If you have no new data, leave this NULL as is by default.
}
}

\details{
The function performs the unconstrained log-contrast regression model as opposed to the log-contrast
regression described in Aitchison (2003), pg. 84-85. The logarithm of the compositional predictor variables
is used (hence no zero values are allowed). The response variable is linked to the log-transformed data
\bold{without} the constraint that the sum of the regression coefficients equals 0. If you want the regression model
with the zum-to-zero contraints see \code{\link{lc.reg}}. Extra predictors variables are allowed as well,
for instance categorical or continuous.
}

\value{
A list including:
\item{be}{
The constrained regression coefficients. Their sum equals 0.
}
\item{covbe}{
If covariance matrix of the constrained regression coefficients.
}
\item{va}{
The estimated regression variance.
}
\item{residuals}{
The vector of residuals.
}
\item{est}{
If the arguments "xnew" and "znew" were given these are the predicted or estimated values, otherwise it is NULL.
}
}

\references{
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
\code{\link{lc.reg}, \link{lcreg.aov}, \link{lc.reg2}, \link{ulc.reg2}, \link{alfa.pcr}, \link{alfa.knn.reg}
}
}

\examples{
y <- iris[, 1]
x <- as.matrix(iris[, 2:4])
x <- x / rowSums(x)
mod1 <- ulc.reg(y, x)
mod2 <- ulc.reg(y, x, z = iris[, 5])
}

