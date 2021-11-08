\name{MLE of distributions defined in the (0, 1) interval}
\alias{beta.est}
\alias{logitnorm.est}
\alias{hsecant01.est}
\alias{simplex.est}
\alias{kumar.est}
\alias{unitweibull.est}
\alias{ibeta.est}
\alias{zilogitnorm.est}
\title{
MLE of distributions defined in the (0, 1) interval
}

\description{
MLE of distributions defined in the (0, 1) interval.
}

\usage{
beta.est(x, tol = 1e-07)
logitnorm.est(x)
hsecant01.est(x, tol = 1e-07)
simplex.est(x, tol = 1e-07)
kumar.est(x, tol = 1e-07)
unitweibull.est(x, tol = 1e-07, maxiters = 100)
ibeta.est(x, tol = 1e-07)
zilogitnorm.est(x)
}

\arguments{
\item{x}{
A numerical vector with proportions, i.e. numbers in (0, 1) (zeros and ones are not allowed).
}
\item{tol}{
The tolerance level up to which the maximisation stops.
}
\item{maxiters}{
The maximum number of iterations the Newton-Raphson algorithm will perform.
}
}

\details{
Maximum likelihood estimation of the parameters of some distributions are performed, some of which use the Newton-Raphson. Some distributions and hence the functions do not accept zeros. "logitnorm.mle" fits the logistic normal, hence no Newton-Raphson is required and the "hypersecant01.mle" and "simplex.est" use the golden ratio search as is it faster than the Newton-Raphson (less computations).
The "zilogitnorm.est" stands for the zero inflated logistic normal distribution.
The "ibeta.est" fits the zero or the one inflated beta distribution.
}

\value{
A list including:
\item{iters}{
The number of iterations required by the Newton-Raphson.
}
\item{loglik}{
The value of the log-likelihood.
}
\item{param}{
The estimated parameters. In the case of "hypersecant01.est" this is called "theta" as there is only one parameter.
}
}

\references{
Kumaraswamy, P. (1980). A generalized probability density function for double-bounded random processes.
Journal of Hydrology. 46(1-2): 79-88.

Jones, M.C. (2009). Kumaraswamy's distribution: A beta-type distribution with some tractability advantages.
Statistical Methodology. 6(1): 70-81.

Zhang, W. & Wei, H. (2008). Maximum likelihood estimation for simplex distribution nonlinear mixed models
via the stochastic approximation algorithm. The Rocky Mountain Journal of Mathematics, 38(5): 1863-1875.

You can also check the relevant wikipedia pages.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{diri.nr2},
}
}

\examples{
x <- rbeta(1000, 1, 4)
beta.est(x)
ibeta.est(x)

x <- runif(1000)
hsecant01.est(x)
logitnorm.est(x)
ibeta.est(x)

x <- rbeta(1000, 2, 5)
x[sample(1:1000, 50)] <- 0
ibeta.est(x)
}
