\name{Log-contrast logistic or Poisson regression with compositional predictor variables}
\alias{lc.glm}
\title{
Log-contrast logistic or Poisson regression with compositional predictor variables
}

\description{
Log-contrast logistic or Poisson regression with compositional predictor variables.
}

\usage{
lc.glm(y, x, z = NULL, model = "logistic", xnew = NULL, znew = NULL)
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
\item{model}{
This can be either "logistic" or "poisson".
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
The function performs the log-contrast logistic or Poisson regression model. The logarithm of the
compositional predictor variables is used (hence no zero values are allowed). The response variable
is linked to the log-transformed data with the constraint that the sum of the regression coefficients
equals 0. If you want the regression without the zum-to-zero contraints see \code{\link{ulc.glm}}.
Extra predictors variables are allowed as well, for instance categorical or continuous.
}

\value{
A list including:
\item{runtime}{
The duration of the algorithm required to minimize the deviance.
}
\item{devi}{
The residual deviance of the logistic or Poisson regression model.
}
\item{be}{
The constrained regression coefficients. Their sum equals 0.
}
\item{est}{
If the arguments "xnew" and znew were given these are the predicted or estimated values, otherwise it is NULL.
}
}

\references{
Aitchison J. (1986). The statistical analysis of compositional data. Chapman & Hall.

Lu J., Shi P., and Li H. (2019). Generalized linear models with linear constraints
for microbiome compositional data. Biometrics, 75(1): 235-244.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{ulc.glm}, \link{lc.reg}, \link{lc.reg2}, \link{alfa.pcr}, \link{glm.pcr}
}
}

\examples{
y <- rbinom(150, 1, 0.5)
x <- rdiri(150, runif(3, 1, 4) )
x <- x / rowSums(x)
mod1 <- lc.glm(y, x)
}

