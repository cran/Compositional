\name{Multivariate kernel density estimation for compositional data}
\alias{comp.kern}
\title{
Multivariate kernel density estimation for compositional data
}

\description{
Multivariate kernel density estimation for compositional data.
}

\usage{
comp.kern(x, type= "alr", h = NULL, thumb = "silverman")
}

\arguments{
\item{x}{
A matrix with Euclidean (continuous) data.
}
\item{type}{
The type of trasformation used, either the additive log-ratio ("alr"), the
isometric log-ratio ("ilr") or the pivot coordinate ("pivot") transformation.
}
\item{h}{
The bandwidh value. It can be a single value, which is turned into a vector and
then into a diagonal matrix, or a vector which is turned into a diagonal matrix.
If it is NULL, then you need to specify the "thumb" argument below.
}
\item{thumb}{
Do you want to use a rule of thumb for the bandwidth parameter? If no, leave
the "h" NULL and put "estim" for maximum likelihood cross-validation, "scott"
or "silverman" for Scott's and Silverman's rules of thumb respectively.
}
}

\details{
The multivariate kernel density estimate is calculated with a (not necssarily
given) bandwidth value.
}

\value{
A vector with the density estimates calculated for every vector.
}

\references{
Arsalane Chouaib Guidoum (2015). Kernel Estimator and Bandwidth Selection for
Density and its Derivatives.

The kedd R package.

M.P. Wand and M.C. Jones (1995). Kernel smoothing, pages 91-92.

B.W. Silverman (1986). Density estimation for statistics and data analysis, pages 76-78.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}
\seealso{
\code{\link{comp.kerncontour}, \link{mkde}
}
}

\examples{
x <- as.matrix(iris[, 1:3])
x <- x / rowSums(x)
f <- comp.kern(x)
}

