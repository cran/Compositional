\name{Bayesian network learning with compositional data}
\alias{compbn}
\title{
Bayesian network learning with compositional data
}

\description{
Bayesian network learning with compositional data.
}


\usage{
compbn(x, type = "fedhc", max_k = 3, alpha = 0.05, robust = FALSE,
ini.stat = NULL, R = NULL, restart = 10, tabu = 10, score = "bic-g",
blacklist = NULL, whitelist = NULL)
}

\arguments{
\item{x}{
A numerical matrix with the compositional data. They can be either the logged compositional or the 
centred log-ratio transformed compositional data. We leave this open to the user.
}
\item{type}{
This can be either "fedhc", "pchc", "mmhc", "fedtabu", "pctabu" or "mmtabu".
}
\item{max_k}{
The maximum conditioning set to use in the conditional indepedence test (see Details). 
Integer, default value is 3
}
\item{alpha}{
The significance level for assessing the p-values.
}
\item{robust}{
Do you want outliers to be removed prior to applying the algorithms? If yes, set this to TRUE to u
tilise the MCD.
}
\item{ini.stat}{
If the initial test statistics (univariate associations) are available, pass them through this parameter.
}
\item{R}{
If the correlation matrix is available, pass it here.
}
\item{restart}{
An integer, the number of random restarts.
}
\item{tabu}{
An integer, the length of the tabu list used in the tabu function.
}
\item{score}{
A character string, the label of the network score to be used in the algorithm. If none is specified,
the default score is the Bayesian Information Criterion. Other available scores are: "bic-g" (default), 
"loglik-g", "aic-g", "bic-g" or "bge".
}
\item{blacklist}{
A data frame with two columns (optionally labeled "from" and "to"), containing a set of arcs not to
be included in the graph.
}
\item{whitelist}{
A data frame with two columns (optionally labeled "from" and "to"), containing a set of arcs to be
included in the graph.
}
}

\details{
The FEDHC algorithm is implemented. The FBED algortihm (Borboudakis and Tsamardinos, 2019), without 
the backward phase, is implemented during the skeleton identification phase. Next, the Hill Climbing 
greedy search or the Tabu search is employed to score the network.

The PC algorithm as proposed by Spirtes et al. (2001) is first implemented followed by a scoring phase, 
such as hill climbing or tabu search. The PCHC was proposed by Tsagris (2021), while the PCTABU algorithm 
is the same but instead of the hill climbing scoring phase, the tabu search is employed.

The MMHC algorithm is implemented without performing the backward elimination during the skeleton 
identification phase. The MMHC as described in Tsamardinos et al. (2006) employs the MMPC algorithm 
during the skeleton construction phase and the Tabu search in the scoring phase. In this package, 
the mmhc function employs the Hill Climbing greedy search in the scoring phase while the mmtabu 
employs the Tabu search.
}

\value{
A list including:
\item{ini}{
A list including the output of the \code{ \link{mmhc.skel} } function.
}
\item{dag}{
A "bn" class output. A list including the outcome of the Hill-Climbing or the Tabu search phase. 
See the package "bnlearn" for more details.
}
\item{scoring}{
The score value.
}
\item{runtime}{
The duration of the algorithm.
}
}

\references{
Tsagris M. (2021). A new scalable Bayesian network learning algorithm with applications to economics. 
Computational Economics, 57(1):341-367.

Tsagris M. (2021). The FEDHC Bayesian network learning algorithm. https://arxiv.org/pdf/2012.00113.pdf.

Borboudakis G. and Tsamardinos I. (2019). Forward-backward selection with early dropping. 
Journal of Machine Learning Research, 20(8): 1-39.

Tsamardinos I., Brown E.L. and Aliferis F.C. (2006). The max-min hill-climbing Bayesian network structure 
learning algorithm. Machine Learning, 65(1): 31-78.

Spirtes P.,  Glymour C. and Scheines R. (2001). Causation, Prediction, and Search. 
The MIT Press, Cambridge, MA, USA, 3nd edition.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}
\seealso{
\code{ \link{acor}, \link{alr}, \link{alfa} }
}

\examples{
# simulate a dataset with continuous data
x <- rdiri( 100, runif(20) )
a <- compbn( log(x) )
}

