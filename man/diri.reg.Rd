\name{Dirichlet regression}
\alias{diri.reg}
\alias{diri.reg2}
\title{
Dirichlet regression
}

\description{
Dirichlet regression.
}

\usage{
diri.reg(y, x, plot = TRUE, xnew = NULL)

diri.reg2(y, x, xnew = NULL)
}

\arguments{
\item{y}{
A matrix with the compositional data (dependent variable). Zero values are not allowed.
}
\item{x}{
The predictor variable(s), they can be either continnuous or categorical or both.
}
\item{plot}{
A boolean variable specifying whether to plot the leverage values of the observations or not. 
This is taken into account only when xnew = NULL.
}
\item{xnew}{
If you have new data use it, otherwise leave it NULL.
}
}

\details{
A Dirichlet distribution is assumed for the regression. This involves numerical optimisation. 
The function "diri.reg2" allows for the covariates to be linked with the precisio parameter 
\eqn{\phi} via the exponential link function \eqn{\phi = e^{x*b}}.
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
\item{phipar}{
The coefficients of the phi parameter if it is linked to the covariates.
}
\item{std.phi}{
The standard errors of the coefficients of the phi parameter is it linked to the covariates.
}
\item{log.phi}{
The logarithm of the precision parameter.
}
\item{std.logphi}{
The standard error of the logarithm of the precision parameter.
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
  \item{lev}{
The leverage values.
}
\item{est}{
For the "diri.reg" this contains the fitted or the predicted values (if xnew is not NULL). 
For the "diri.reg2" if xnew is NULL, this
is also NULL.
}
}

\references{
Maier, Marco J. (2014) DirichletReg: Dirichlet Regression for Compositional Data in R.
Research Report Series/Department of Statistics and Mathematics, 125. WU Vienna University of 
Economics and Business, Vienna.
http://epub.wu.ac.at/4077/1/Report125.pdf

Gueorguieva, Ralitza, Robert Rosenheck, and Daniel Zelterman (2008). Dirichlet component 
regression and its applications to psychiatric data. Computational statistics & data analysis 
52(12): 5344-5355.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou 
<gioathineou@gmail.com>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{js.compreg}, \link{kl.compreg}, \link{ols.compreg}, \link{comp.reg}, \link{alfa.reg}
}
}

\examples{
x <- as.vector(iris[, 4])
y <- as.matrix(iris[, 1:3])
y <- y / rowSums(y)
mod1 <- diri.reg(y, x)
mod2 <-diri.reg2(y, x)
mod3 <- comp.reg(y, x)
}

\keyword{ Dirichlet regression }
\keyword{ multivariate regression }
