\name{Projection pursuit regression with compositional predictor variables}
\alias{pprcomp}
\title{
Projection pursuit regression with compositional predictor variables
}

\description{
Projection pursuit regression with compositional predictor variables.
}

\usage{
pprcomp(y, x, nterms = 3, type = "alr", xnew = NULL)
}

\arguments{
\item{y}{
A numerical vector with the continuous variable.
}
\item{x}{
A matrix with the compositional data. No zero values are allowed.
}
\item{nterms}{
The number of terms to include in the final model.
}
\item{type}{
Either "alr" or "log" corresponding to the additive log-ratio transformation or the simple
logarithm applied to the compositional data.
}
\item{xnew}{
If you have new data use it, otherwise leave it NULL.
}
}

\details{
This is the standard projection pursuit. See the built-in function "ppr" for more details.
When the data are transformed with the additive log-ratio transformation this is close
in spirit to the log-contrast regression.
}

\value{
A list includign:
\item{runtime}{
The runtime of the regression.
}
\item{mod}{
The produced model as returned by the function "ppr".
}
\item{est}{
The fitted values of xnew if xnew is not NULL.
}
}

\references{
Friedman, J. H. and Stuetzle, W. (1981). Projection pursuit regression. Journal of the American
Statistical Association, 76, 817-823. doi: 10.2307/2287576.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{pprcomp.tune}, \link{ice.pprcomp}, \link{alfa.pcr}, \link{lc.reg}, \link{comp.ppr}
}
}

\examples{
x <- as.matrix( iris[, 2:4] )
x <- x/ rowSums(x)
y <- iris[, 1]
pprcomp(y, x)
}


