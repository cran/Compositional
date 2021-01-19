\name{Unit-Weibull regression models for proportions}
\alias{unitweib.reg}
\title{
Unit-Weibull regression models for proportions
}

\description{
Unit-Weibull regression models for proportions.
}

\usage{
unitweib.reg(y, x, tau = 0.5)
}

\arguments{
\item{y}{
A numerical vector proportions. 0s and 1s are allowed.
}
\item{x}{
A matrix or a data frame with the predictor variables.
}
\item{tau}{
The quantile to be used for estimation. The default value is 0.5 yielding the median.
}
}

\details{
See the reference paper.
}

\value{
A list including:
\item{loglik}{
The loglikelihood of the regression model.
}
\item{info}{
A matrix with all estimated parameters, their standard error, their Wald-statistic and its associated p-value.
}
}

\references{
Mazucheli J., Menezes A. F. B., Fernandes L. B., de Oliveira R. P. and Ghitany M. E. (2020).
The unit-Weibull distribution as an alternative to the Kumaraswamy distribution for the modeling of quantiles conditional on covariates.
Journal of Applied Statistics, 47(6): 954--974.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{propreg}, \link{beta.reg}
}
}

\examples{
y <- exp( - rweibull(100, 1, 1) )
x <- matrix( rnorm(100 * 2), ncol = 2 )
a <- unitweib.reg(y, x)
}

