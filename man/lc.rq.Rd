\name{Log-contrast quantile regression with compositional predictor variables}
\alias{lc.rq}
\title{
Log-contrast quantile regression with compositional predictor variables
}

\description{
Log-contrast quantile regression with compositional predictor variables.
}

\usage{
lc.rq(y, x, z = NULL, tau, xnew = NULL, znew = NULL)
}

\arguments{
\item{y}{
A numerical vector containing the response variable values.
}
\item{x}{
A matrix with the predictor variables, the compositional data.
No zero values are allowed.
}
\item{z}{
A matrix, data.frame, factor or a vector with some other covariate(s).
}
\item{tau}{
The quantile to be estimated, a number between 0 and 1.
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
The function performs the quantile regression model. The logarithm of the compositional
predictor variables is used (hence no zero values are allowed). The response variable is
linked to the log-transformed data with the constraint that the sum of the regression
coefficients equals 0. If you want the regression without the zum-to-zero contraints see \code{\link{ulc.rq}}. 
Extra predictor variables are allowed as well, for instance categorical
or continuous.
}

\value{
A list including:
\item{mod}{
The object as returned by the function quantreg::rq(). This is useful for hypothesis
testing purposes.
}
\item{be}{
The constrained regression coefficients. Their sum (excluding the constant) equals 0.
}
\item{est}{
If the arguments "xnew" and znew were given these are the predicted or estimated values,
otherwise it is NULL.
}
}

\references{
Aitchison J. (1986). The statistical analysis of compositional data. Chapman & Hall.

Koenker R. W. and Bassett G. W. (1978). Regression Quantiles, Econometrica, 46(1): 33--50.

Koenker R. W. and d'Orey V. (1987). Algorithm AS 229: Computing Regression Quantiles.
Applied Statistics, 36(3): 383--393.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{lc.rq2}, \link{ulc.rq}
}
}

\examples{
y <- rnorm(150)
x <- rdiri(150, runif(3, 1, 4) )
mod1 <- lc.rq(y, x)
}

