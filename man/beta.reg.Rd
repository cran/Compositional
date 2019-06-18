\name{Beta regression}
\alias{beta.reg}
\title{
Beta regression
}

\description{
Beta regression.
}

\usage{
beta.reg(y, x, xnew = NULL)
}

\arguments{
\item{y}{
The response variable. It must be a numerical vector with proportions excluding 0 and 1.
}
\item{x}{
The indendent variable(s). It can be a vector, a matrix or a dataframe with continuous only variables,
a data frame with mixed or only categorical variables.
}
\item{xnew}{
If you have new values for the predictor variables (dataset) whose response values you want to predict insert them here.
}
}

\details{
Beta regression is fitted.
}

\value{
A list including:
\item{phi}{
The estimated precision parameter.
}
\item{info}{
A matrix with the estimated regression parameters, their standard errors, Wald statistics and associated p-values.
}
\item{loglik}{
The log-likelihood of the regression model.
}
\item{est}{
The estimated values if xnew is not NULL.
}
}

\references{
Ferrari S.L.P. and Cribari-Neto F. (2004). Beta Regression for Modelling Rates and Proportions.
Journal of Applied Statistics, 31(7): 799-815.
}

\author{
 Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{beta.est}, \link{prop.reg}, \link{diri.reg} }
}

\examples{
y <- rbeta(300, 3, 5)
x <- matrix( rnorm(300 * 2), ncol = 2)
beta.reg(y, x)
}

