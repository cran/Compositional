\name{Zero adjusted Dirichlet regression}
\alias{zadr}
\alias{mixreg}
\title{
Zero adjusted Dirichlet regression
}

\description{
Zero adjusted Dirichlet regression.
}

\usage{
zadr(y, x, xnew = NULL, tol = 1e-05)

mixreg(param, z)
}

\arguments{
\item{y}{
A matrix with the compositional data (dependent variable). The number of observations with at least one zero value should not be more
than the columns of the predictor variables. Otherwise, the initial values will not be calculated.
}
\item{x}{
The predictor variable(s), they can be either continnuous or categorical or both.
}
\item{xnew}{
If you have new data use it, otherwise leave it NULL.
}
\item{tol}{
A tolerance level to terminate the maximisation process.
}
\item{param}{
Some arguments passed on to the mixreg helper function.
}
\item{z}{
Some arguments passed on to the mixreg helper function.
}
}

\details{
A zero adjusted Dirichlet regression is being fittd. The likelihood conists of two components. The contributions of the non zero
compositional values and the contributions of the compositional vectors with at least one zero value. The second component may
have many different sub-categories, one for each pattern of zeros.
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
The precision parameter. If covariates are linked with it (function "diri.reg2"), this will
be a vector.
}
\item{be}{
The beta coefficients.
}
\item{seb}{
The standard error of the beta coefficients.
}
\item{sigma}{
Th covariance matrix of the regression parameters (for the mean vector and the phi parameter) in
the function "diri.reg2".
}
\item{est}{
The fitted or the predicted values (if xnew is not NULL).
}
}

\references{
Tsagris M. and Stewart C. (2018). A Dirichlet Regression Model for Compositional Data with Zeros.
Accepted at the Lobachevskii Journal of Mathematics.

Preprint available from https://arxiv.org/pdf/1410.5011.pdf
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{diri.reg}, \link{kl.compreg}, \link{ols.compreg}, \link{alfa.reg}
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

