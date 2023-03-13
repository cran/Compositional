\name{Zero adjusted Dirichlet regression}
\alias{zadr}
\alias{zadr2}
\title{
Zero adjusted Dirichlet regression
}

\description{
Zero adjusted Dirichlet regression.
}

\usage{
zadr(y, x, con = TRUE, B = 1, ncores = 2, xnew = NULL)
zadr2(y, x, con = TRUE, B = 1, ncores = 2, xnew = NULL)
}

\arguments{
\item{y}{
A matrix with the compositional data (dependent variable). The number of observations
(vectors) with no zero values should be more than the columns of the predictor variables.
Otherwise, the initial values will not be calculated.
}
\item{x}{
The predictor variable(s), they can be either continnuous or categorical or both.
}
\item{con}{
If this is TRUE (default) then the constant term is estimated, otherwise the model includes no constant term.
}
\item{B}{
If B is greater than 1 bootstrap estimates of the standard error are returned.
If you set this greater than 1, then you must define the number of clusters in
order to run in parallel.
}
\item{ncores}{
The number of cores to use when B>1. This is to be used for the
case of bootstrap. If B = 1, this is not taken into consideration.
If this does not work then you might need to load the doParallel yourselves.
}
\item{xnew}{
If you have new data use it, otherwise leave it NULL.
}
}

\details{
A zero adjusted Dirichlet regression is being fittd. The likelihood conists of two components.
The contributions of the non zero compositional values and the contributions of the compositional
vectors with at least one zero value. The second component may have many different sub-categories,
one for each pattern of zeros. The function "zadr2()" links the covariates to the alpha parameters
of the Dirichlet distribution, i.e. it uses the classical parametrization of the distribution.
This means, that there is a set of regression parameters for each component.
}

\value{
A list including:
\item{runtime}{
The time required by the regression.
}
\item{loglik}{
The value of the log-likelihood.
}
\item{phi}{
The precision parameter.
}
\item{be}{
The beta coefficients.
}
\item{seb}{
The standard error of the beta coefficients.
}
\item{sigma}{
Th covariance matrix of the regression parameters (for the mean vector and the phi parameter).
}
\item{est}{
The fitted or the predicted values (if xnew is not NULL).
}
}

\references{
Tsagris M. and Stewart C. (2018). A Dirichlet regression model for compositional data with zeros.
Lobachevskii Journal of Mathematics,39(3): 398--412.

Preprint available from https://arxiv.org/pdf/1410.5011.pdf
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{zad.est}, \link{diri.reg}, \link{kl.compreg}, \link{ols.compreg}, \link{alfa.reg}
}
}

\examples{
x <- as.vector(iris[, 4])
y <- as.matrix(iris[, 1:3])
y <- y / rowSums(y)
mod1 <- diri.reg(y, x)
y[sample(1:450, 15) ] <- 0
mod2 <- zadr(y, x)
}

