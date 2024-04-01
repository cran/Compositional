\name{Ternary diagram with confidence region for the matrix of coefficients of the SCLS or the TFLR model}
\alias{ternary.coefcr}

\title{
Ternary diagram with confidence region for the matrix of coefficients of the SCLS or the TFLR model
}

\description{
Ternary diagram with confidence region for the matrix of coefficients of the SCLS or the TFLR model.
}

\usage{
ternary.coefcr(y, x, type = "scls", conf = 0.95, R = 1000, dg = FALSE, hg = FALSE)
}

\arguments{
\item{y}{
A matrix with the response compositional data.
}
\item{x}{
A matrix with the predictor compositional data.
}
\item{type}{
The type of model to use, "scls" or "tflr". Depending on the model selected, the function will construct the confidence regions of the estimated matrix of coefficients of that model.
}
\item{conf}{
The confidence level, by default this is set to 0.95.
}
\item{R}{
Number of bootstrap replicates to run.
}
\item{dg}{
Do you want diagonal grid lines to appear? If yes, set this TRUE.
}
\item{hg}{
Do you want horizontal grid lines to appear? If yes, set this TRUE.
}
}

\details{
This function runs the SCLS or the TFLR model and constructs confidence regions for the estimated matrix of regression coefficients using non-parametric bootstrap.
}

\value{
A ternary plot of the estimated matrix of coefficients of the SCLS or of the TFLR model, and their associated confidence regions.
}

\references{
Fiksel J., Zeger S. and Datta A. (2022). A transformation-free linear regression for compositional outcomes and predictors. Biometrics, 78(3): 974--987.

Tsagris. M. (2024). Constrained least squares simplicial-simplicial regression.
https://arxiv.org/pdf/2403.19835.pdf
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{ternary}, \link{scls}, \link{tflr}, \link{ternary.mcr}
}
}

\examples{
y <- rdiri(50, runif(3))
x <- rdiri(50, runif(4))
ternary.coefcr(y, x, R = 500, dg = TRUE, hg = TRUE)
}
