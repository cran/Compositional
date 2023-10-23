\name{ANOVA for the log-contrast GLM versus the uncostrained GLM}
\alias{lcglm.aov}
\title{
ANOVA for the log-contrast GLM versus the uncostrained GLM
}

\description{
ANOVA for the log-contrast GLM versus the uncostrained GLM.
}

\usage{
lcglm.aov(mod0, mod1)
}

\arguments{
\item{mod0}{
The log-contrast GLM. The object returned by \code{\link{lc.glm}}.
}
\item{mod1}{
The unconstrained GLM. The object returned by \code{\link{ulc.glm}}.
}
}

\details{
A chi-square test is performed to test the zero-to-sum constraints of the regression coefficients.
}

\value{
A vector with two values, the chi-square test statistic and its associated p-value.
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
\code{\link{lc.glm}, \link{ulc.glm}
}
}

\examples{
y <- rbinom(150, 1, 0.5)
x <- as.matrix(iris[, 2:4])
x <- x / rowSums(x)
mod0 <- lc.glm(y, x)
mod1 <- ulc.glm(y, x)
lcglm.aov(mod0, mod1)
}

