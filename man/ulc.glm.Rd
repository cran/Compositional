\name{Unconstrained GLMs with compositional predictor variables}
\alias{ulc.glm}
\title{
Unconstrained GLMs with compositional predictor variables
}

\description{
Unconstrained GLMs with compositional predictor variables.
}

\usage{
ulc.glm(y, x, z = NULL, model = "logistic", xnew = NULL, znew = NULL)
}

\arguments{
\item{y}{
A numerical vector containing the response variable values. This is either a binary variable or a vector with counts.
}
\item{x}{
A matrix with the predictor variables, the compositional data. No zero values are allowed.
}
\item{z}{
A matrix, data.frame, factor or a vector with some other covariate(s).
}
\item{model}{
For the ulc.glm(), this can be either "logistic" or "poisson".
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
The function performs the unconstrained log-contrast logistic or Poisson regression model. The logarithm of the
compositional predictor variables is used (hence no zero values are allowed). The response variable
is linked to the log-transformed data \bold{without} the constraint that the sum of the regression coefficients
equals 0. If you want the regression without the zum-to-zero contraints see \code{\link{lc.glm}}.
Extra predictors variables are allowed as well, for instance categorical or continuous.
}

\value{
A list including:
\item{devi}{
The residual deviance of the logistic or Poisson regression model.
}
\item{be}{
The unconstrained regression coefficients. Their sum does not equal 0.
}
\item{est}{
If the arguments "xnew" and znew were given these are the predicted or estimated values, otherwise it is NULL.
}
}

\references{
Aitchison J. (1986). The statistical analysis of compositional data. Chapman & Hall.

Lu J., Shi P., and Li H. (2019). Generalized linear models with linear constraints
for microbiome compositional data. Biometrics, 75(1): 235--244.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{lc.glm}, \link{lc.glm2}, \link{ulc.glm2},  \link{lcglm.aov}
}
}

\examples{
y <- rbinom(150, 1, 0.5)
x <- rdiri(150, runif(3, 1,3))
mod <- ulc.glm(y, x)
}

