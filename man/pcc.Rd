\name{Proportionality correlation coefficient matrix}
\alias{pcc}
\title{
Proportionality correlation coefficient matrix
}

\description{
Proportionality correlation coefficient matrix.
}

\usage{
pcc(x)
}

\arguments{
\item{x}{
A numerical matrix with the compositional data. Zeros are not allowed as the logarithm is applied.
}
}

\details{
The function returns the proportionality correlation coefficient matrix. See Lovell et al. (2015) for more information.
}

\value{
A matrix with the alr transformed data (if alr is used) or with the compositional data (if the alrinv is used).
}

\references{
Zheng, B. (2000). Summarizing the goodness of fit of generalized linear models for longitudinal data.
Statistics in medicine, 19(10), 1265-1275.

Lovell D., Pawlowsky-Glahn V., Egozcue J. J., Marguerat S. and Bahler, J. (2015). Proportionality: a valid alternative
to correlation for relative data. PLoS Computational Biology, 11(3), e1004075.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{acor}, \ link{alr}
}
}

\examples{
x <- Compositional::rdiri(100, runif(4) )
a <- Compositional::pcc(x)
}

