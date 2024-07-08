\name{MLE of the zero adjusted Dirichlet distribution}
\alias{zad.est}
\title{
MLE of the zero adjusted Dirichlet distribution
}

\description{
MLE of the zero adjusted Dirichlet distribution.
}

\usage{
zad.est(y)
}

\arguments{
\item{y}{
A matrix with the compositional data.
}
}

\details{
A zero adjusted Dirichlet distribution is being fitted and its parameters are estimated.
}

\value{
A list including:
\item{loglik}{
The value of the log-likelihood.
}
\item{phi}{
The precision parameter. If covariates are linked with it (function "diri.reg2"), this will
be a vector.
}
\item{mu}{
The mean vector of the distribution.
}
\item{runtime}{
The time required by the regression.
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
\code{\link{zadr}, \link{diri.nr}, \link{zilogitnorm.est}, \link{zeroreplace}
}
}

\examples{
y <- as.matrix(iris[, 1:3])
y <- y / rowSums(y)
mod1 <- diri.nr(y)
y[sample(1:450, 15) ] <- 0
mod2 <- zad.est(y)
}

