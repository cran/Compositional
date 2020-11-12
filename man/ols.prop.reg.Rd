\name{Distance based regression models for proportions}
\alias{ols.prop.reg}
\alias{js.prop.reg}
\alias{helling.prop.reg}
\title{
Distance based regression models for proportions
}

\description{
Distance based regression models for proportions.
}

\usage{
ols.prop.reg(y, x, cov = FALSE, tol = 1e-07, maxiters = 100)
js.prop.reg(y, x, tol = 1e-07, maxiters = 100)
helling.prop.reg(y, x, tol = 1e-07, maxiters = 100)
}

\arguments{
\item{y}{
A numerical vector proportions. 0s and 1s are allowed.
}
\item{x}{
A matrix or a data frame with the predictor variables.
}
\item{cov}{
Should the covariance matrix be returned? TRUE or FALSE.
}
\item{tol}{
The tolerance value to terminate the Newton-Raphson algorithm. This is set to \eqn{10^{-9}} by default.
}
\item{maxiters}{
The maximum number of iterations before the Newton-Raphson is terminated automatically.
}
}

\details{
We are using the Newton-Raphson, but unlike R's built-in function "glm" we do no checks and no extra calculations, or whatever. Simply the model. The functions accept binary responses as well (0 or 1).
}

\value{
A list including:
\item{sse}{
The sum of squres of errors for the "ols.prop.reg" function.
}
\item{be}{
The estimated regression coefficients.
}
\item{seb}{
The standard error of the regression coefficients in "ols.prop.reg" if "cov" is TRUE.
}
\item{covb}{
The covariance matrix of the regression coefficients in "ols.prop.reg" if "cov" is TRUE.
}
\item{js}{
The Jensen-Shannon divergence between the true and the obseervd proportions in "js.prop.reg".
}
\item{H}{
The Hellinger distance between the true and the obseervd proportions in "helling.prop.reg".
}
\item{iters}{
The number of iterations required by the Newton-Raphson.
}
}

\references{
Papke L. E. & Wooldridge J. (1996). Econometric methods for fractional response variables with
an application to 401(K) plan participation rates. Journal of Applied Econometrics, 11(6): 619--632.

McCullagh, Peter, and John A. Nelder. Generalized linear models. CRC press, USA, 2nd edition, 1989.
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
y <- rbeta(100, 1, 4)
x <- matrix(rnorm(100 * 2), ncol = 2)
a1 <- ols.prop.reg(y, x)
a2 <- js.prop.reg(y, x)
a3 <- helling.prop.reg(y, x)
}

