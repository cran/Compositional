\name{Log-contrast regression with multiple compositional predictors}
\alias{lc.reg2}
\title{
Log-contrast regression with multiple compositional predictors
}

\description{
Log-contrast regression with multiple compositional predictors.
}

\usage{
lc.reg2(y, x, z = NULL, xnew = NULL, znew = NULL)
}

\arguments{
\item{y}{
A numerical vector containing the response variable values. This must be a continuous variable.
}
\item{x}{
A list with multiple matrices with the predictor variables, the compositional data. No zero values are allowed.
}
\item{z}{
A matrix, data.frame, factor or a vector with some other covariate(s).
}
\item{xnew}{
A matrix containing a list with multiple matrices with compositional data whose response is to be predicted.
If you have no new data, leave this NULL as is by default.
}
\item{znew}{
A matrix, data.frame, factor or a vector with the values of some other covariate(s).
If you have no new data, leave this NULL as is by default.
}
}

\details{
The function performs the log-contrast regression model as described in Aitchison (2003), pg. 84-85.
The logarithm of the compositional predictor variables is used (hence no zero values are allowed).
The response variable is linked to the log-transformed data with the constraint that the sum of the
regression coefficients equals 0. Hence, we apply constrained least squares, which has a closed form
solution. The constrained least squares is described in Chapter 8.2 of Hansen (2019). The idea is to
minimise the sum of squares of the residuals under the constraint \eqn{R^T \beta = c}, where \eqn{c=0}
in our case. If you want the regression without the zum-to-zero contraints see \code{\link{ulc.reg2}}.
Extra predictors variables are allowed as well, for instance categorical or continuous. The difference
with \code{\link{lc.reg}} is that instead of one, there are multiple compositions treated as predictor
variables.
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

Hansen, B. E. (2022). Econometrics. Princeton University Press.

Xiaokang Liu, Xiaomei Cong, Gen Li, Kendra Maas and Kun Chen (2020). Multivariate Log-Contrast Regression
with Sub-Compositional Predictors: Testing the Association Between Preterm Infants' Gut Microbiome and
Neurobehavioral Outcome.
\url{https://arxiv.org/pdf/2006.00487.pdf}.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{ulc.reg2}, \link{lc.reg}, \link{ulc.reg}, \link{lcreg.aov}, \link{alfa.pcr}, \link{alfa.knn.reg}
}
}

\examples{
y <- iris[, 1]
x <- list()
x1 <- as.matrix(iris[, 2:4])
x1 <- x1 / rowSums(x1)
x[[ 1 ]] <- x1
x[[ 2 ]] <- rdiri(150, runif(4) )
x[[ 3 ]] <- rdiri(150, runif(5) )
mod <- lc.reg2(y, x)
}

