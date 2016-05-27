\name{diri.reg}
\alias{diri.reg}
\title{
Dirichlet regression
}

\description{
Dirichlet regression.
}

\usage{
diri.reg(y, x, plot = TRUE, xnew = NULL)
}

\arguments{
  \item{y}{
  A matrix with the compositional data (dependent variable). Zero values are not allowed.
}
  \item{x}{
  The predictor variable(s), they can be either continnuous or categorical or both.
}
  \item{plot}{
  A boolean variable specifying whether to plot the leverage values of the observations or not. This is taken into account only when xnew = NULL.
}
  \item{xnew}{
  If you have new data use it, otherwise leave it NULL.
}
}

\details{
A Dirichlet distribution is assumed for the regression. This involves numerical optimisation.
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
  \item{log.phi}{
The logarithm of the precision parameter.
}
  \item{std.logphi}{
The standard error of the logarithm of the precision parameter.
}
  \item{beta}{
The beta coefficients.
}
  \item{seb}{
The standard error of the beta coefficients.
}
  \item{lev}{
The leverage values.
}
  \item{est}{
The fitted or the predicted values (if xnew is not NULL).
}
}

\references{
Maier, Marco J. (2014) DirichletReg: Dirichlet Regression for Compositional Data in R.
Research Report Series/Department of Statistics and Mathematics, 125. WU Vienna University of Economics and Business, Vienna.
http://epub.wu.ac.at/4077/1/Report125.pdf

Gueorguieva, Ralitza, Robert Rosenheck, and Daniel Zelterman (2008). Dirichlet component regression and its applications to psychiatric data.
Computational statistics & data analysis 52(12): 5344-5355.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{js.compreg}, \link{kl.compreg}, \link{ols.compreg}, \link{comp.reg}, \link{alfa.reg}
}
}

\examples{
x <- iris[, 4]
y <- iris[, 1:3]
mod1 <- diri.reg(y, x)
mod2 <- comp.reg(y, x)
}

\keyword{ Dirichlet regression }
\keyword{ multivariate regression }
