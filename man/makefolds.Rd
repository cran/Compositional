\name{Generate random folds for cross-validation}
\alias{makefolds}
\title{
Generate random folds for cross-validation
}

\description{
Random folds for use in a cross-validation are generated. There is the option for stratified splitting as well.
}

\usage{
makefolds(ina, nfolds = 10, stratified = TRUE, seed = NULL)
}

\arguments{
\item{ina}{
A variable indicating the groupings.
}
\item{nfolds}{
The number of folds to produce.
}
\item{stratified}{
A boolean variable specifying whether stratified random (TRUE) or simple random (FALSE) sampling is to be used when producing the folds.
}
\item{seed}{
You can specify your own seed number here or leave it NULL.
}
}

\details{
I was inspired by the command in the package \bold{TunePareto} in order to do the stratified version.
}


\value{
A list with nfolds elements where each elements is a fold containing the indices of the data.
}

%\references{
%}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

\seealso{
\code{\link{compknn.tune} }
}

\examples{
a <- makefolds(iris[, 5], nfolds = 5, stratified = TRUE)
table(iris[a[[1]], 5])  ## 10 values from each group
}


