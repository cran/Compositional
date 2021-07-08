\name{ANOVA for the log-contrast regression versus the uncostrained log-contrast regression}
\alias{lcreg.aov}
\title{
ANOVA for the log-contrast regression versus the uncostrained log-contrast regression
}

\description{
ANOVA for the log-contrast regression versus the uncostrained log-contrast regression.
}

\usage{
lcreg.aov(mod0, mod1)
}

\arguments{
\item{mod0}{
The log-contrast regression model. The object returned by \code{\link{lc.reg}}.
}
\item{mod1}{
The unconstrained log-contrast regression model. The object returned by \code{\link{ulc.reg}}.
}
}

\details{
An F-test is performed to test the zero-to-sum constraints of the regression coefficients.
}

\value{
A vector with two vaues, the F-statistic and its associated p-value.
}

%\references{
%Aitchison J. (1986). The statistical analysis of compositional data. Chapman & Hall.
%}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{lc.reg}, \link{ulc.reg}, \link{alfa.pcr}, \link{alfa.knn.reg}
}
}

\examples{
y <- iris[, 1]
x <- as.matrix(iris[, 2:4])
x <- x / rowSums(x)
mod0 <- lc.reg(y, x)
mod1 <- ulc.reg(y, x)
lcreg.aov(mod0, mod1)
}

