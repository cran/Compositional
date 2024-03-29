\name{Multivariate kernel density estimation}
\alias{mkde}
\title{
Multivariate kernel density estimation
}

\description{
Multivariate kernel density estimation.
}

\usage{
mkde(x, h = NULL, thumb = "silverman")
}

\arguments{
\item{x}{
A matrix with Euclidean (continuous) data.
}
\item{h}{
The bandwidh value. It can be a single value, which is turned into a vector and
then into a diagonal matrix, or a vector which is turned into a diagonal matrix.
If you put this NULL then you need to specify the "thumb" argument below.
}
\item{thumb}{
Do you want to use a rule of thumb for the bandwidth parameter? If no, set h
equal to NULL and put "estim" for maximum likelihood cross-validation, "scott"
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
Arsalane Chouaib Guidoum (2015). Kernel Estimator and Bandwidth Selection for Density and its Derivatives.
The kedd R package.

M.P. Wand and M.C. Jones (1995). Kernel smoothing, pages 91-92.

B.W. Silverman (1986). Density estimation for statistics and data analysis, pages 76-78.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} and Giorgos Athineou <gioathineou@gmail.com>.
}

%\note{
%%  ~~further notes~~
%}
\seealso{
\code{\link{mkde.tune}, \link{comp.kerncontour}
}
}

\examples{
mkde( as.matrix(iris[, 1:4]), thumb = "scott" )
mkde( as.matrix(iris[, 1:4]), thumb = "silverman" )
}

\keyword{ kernel density estimate }
\keyword{ multivariate kernel }
