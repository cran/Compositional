\name{Exponential empirical likelihood hypothesis testing for two mean vectors}
\alias{eel.test2}
\title{
Exponential empirical likelihood hypothesis testing for two mean vectors
}

\description{
Exponential empirical likelihood hypothesis testing for two mean vectors.
}

\usage{
eel.test2(y1, y2, tol = 1e-07, R = 0, graph = FALSE)
}

\arguments{
\item{y1}{
A matrix containing the Euclidean data of the first group.
}
\item{y2}{
A matrix containing the Euclidean data of the second group.
}
\item{tol}{
The tolerance level used to terminate the Newton-Raphson algorithm.
}
\item{R}{
If R is 0, the classical chi-square distribution is used, if R = 1,
the corrected chi-square distribution (James, 1954) is used and if R = 2,
the modified F distribution (Krishnamoorthy and Yanping, 2006) is used.
If R is greater than 3 bootstrap calibration is performed.
}
\item{graph}{
A boolean variable which is taken into consideration only when bootstrap
calibration is performed. IF TRUE the histogram of the bootstrap test
statistic values is plotted.
}
}

\details{
Exponential empirical likelihood or exponential tilting was first introduced by Efron (1981) as a way to perform a "tilted" version of the bootstrap for the one sample mean hypothesis testing. Similarly to the empirical likelihood, positive weights \eqn{p_i}, which sum to one, are allocated to the observations, such that the weighted sample mean \eqn{{\bf \bar{x}}} is equal to some population mean \eqn{\pmb{\mu}}, under the \eqn{H_0}. Under \eqn{H_1} the weights are equal to \eqn{\frac{1}{n}}, where \eqn{n} is the sample size. Following Efron (1981), the choice of \eqn{p_is} will minimize the Kullback-Leibler distance from \eqn{H_0} to \eqn{H_1}
\deqn{
D\left(L_0,L_1\right)=\sum_{i=1}^np_i\log\left(np_i\right),
}
subject to the constraint \eqn{\sum_{i=1}^np_i{\bf x}_i=\pmb{\mu}}. The probabilities take the form
\deqn{
p_i=\frac{e^{\pmb{\lambda}^T{\bf x}_i}}{\sum_{j=1}^ne^{\pmb{\lambda}^T{\bf x}_j}}
}
and the constraint becomes
\deqn{
\frac{\sum_{i=1}^ne^{\pmb{\lambda}^T{\bf x}_i}\left({\bf x}_i-\pmb{\mu}\right)}{\sum_{j=1}^ne^{\pmb{\lambda}^T{\bf x}_j}}=0 \Rightarrow \frac{\sum_{i=1}^n{\bf x}_ie^{\pmb{\lambda}^T{\bf x}_i}}{\sum_{j=1}^ne^{\pmb{\lambda}^T{\bf x}_j}}-\pmb{\mu}=0.
}
Similarly to empirical likelihood a numerical search over \eqn{\pmb{\lambda}} is required.

We can derive the asymptotic form of the test statistic in the two sample means case but in a simpler form, generalizing the approach of Jing and Robinson (1997) to the multivariate case as follows. The three constraints are
\deqn{
{\begin{array}{ccc}
\left(\sum_{j=1}^{n_1}e^{\pmb {\lambda}_1^T{\bf x}_j}\right)^{-1}\left(\sum_{i=1}^{n_1}{\bf x}_ie^{\pmb{\lambda}_1^T
{\bf x}_i}\right) -\pmb{\mu} & = & {\bf 0} \\
\left(\sum_{j=1}^{n_2}e^{\pmb {\lambda}_2^T{\bf y}_j}\right)^{-1}\left(\sum_{i=1}^{n_2}{\bf y}_ie^{\pmb{\lambda}_2^T
{\bf y}_i}\right) -\pmb{\mu} & = & {\bf 0} \\
n_1\pmb{\lambda}_1+n_2\pmb{\lambda}_2 & = & {\bf 0}.
\end{array}}
}
Similarly to EL the sum of a linear combination of the \eqn{\pmb{\lambda}s} is set to zero. We can equate the first two constraints of
\deqn{
\left(\sum_{j=1}^{n_1}e^{\pmb {\lambda}_1^T{\bf x}_j}\right)^{-1}\left(\sum_{i=1}^{n_1}{\bf x}_ie^{\pmb{\lambda}_1^T
{\bf x}_i}\right)=
\left(\sum_{j=1}^{n_2}e^{\pmb {\lambda}_2^T{\bf y}_j}\right)^{-1}\left(\sum_{i=1}^{n_2}{\bf y}_ie^{\pmb{\lambda}_2^T
{\bf y}_i}\right).
}
Also, we can write the third constraint of as \eqn{\pmb{\lambda}_2=-\frac{n_1}{n_2}\pmb{\lambda}_1} and thus rewrite the first two constraints as
\deqn{
\left(\sum_{j=1}^{n_1}e^{\pmb{\lambda}^T{\bf x}_j}\right)^{-1}\left(\sum_{i=1}^{n_1}{\bf x}_ie^{\pmb{\lambda}^T
{\bf x}_i}\right) =
\left(\sum_{j=1}^{n_2}e^{-\frac{n_1}{n_2}\pmb{\lambda}^T{\bf y}_j}\right)^{-1}\left(\sum_{i=1}^{n_2}{\bf y}_ie^{-\frac{n_1}{n_2}\pmb{\lambda}^T
{\bf y}_i}\right).
}

This trick allows us to avoid the estimation of the common mean. It is not possible though to do this in the empirical likelihood method. Instead of minimisation of the sum of the one-sample test statistics from the common mean, we can define the probabilities by searching for the \eqn{\pmb{\lambda}} which makes the last equation hold true. The third constraint of is a convenient constraint, but Jing and Robinson (1997) mention that even though as a constraint is simple it does not lead to second-order accurate confidence intervals unless the two sample sizes are equal. Asymptotically, the test statistic follows a \eqn{\chi^2_d} under the null hypothesis.
}

\value{
A list including:
\item{test}{
The empirical likelihood test statistic value.
}
\item{modif.test}{
The modified test statistic, either via the chi-square or the F distribution.
}
\item{dof}{
The degrees of freedom of the chi-square or the F distribution.
}
\item{pvalue}{
The asymptotic or the bootstrap p-value.
}
\item{mu}{
The estimated common mean vector.
}
\item{runtime}{
The runtime of the bootstrap calibration.
}
}

\references{
Efron B. (1981) Nonparametric standard errors and confidence intervals. Canadian Journal of
Statistics, 9(2): 139--158.

Jing B.Y. and Wood A.T.A. (1996). Exponential empirical likelihood is
not Bartlett correctable. Annals of Statistics, 24(1): 365--369.

Jing B.Y. and Robinson J. (1997). Two-sample nonparametric tilting method. Australian Journal of
Statistics, 39(1): 25--34.

Owen A.B. (2001). Empirical likelihood. Chapman and Hall/CRC Press.

Preston S.P. and Wood A.T.A. (2010). Two-Sample Bootstrap Hypothesis Tests
for Three-Dimensional Labelled Landmark Data. Scandinavian Journal of
Statistics 37(4): 568--587.

Tsagris M., Preston S. and Wood A.T.A. (2017). Nonparametric hypothesis
testing for equality of means on the simplex.
Journal of Statistical Computation and Simulation, 87(2): 406--422.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{el.test2}, \link{maovjames}, \link{hotel2T2},
\link{james}
}
}

\examples{
y1 = as.matrix(iris[1:25, 1:4])
y2 = as.matrix(iris[26:50, 1:4])
eel.test2(y1, y2)
eel.test2(y1, y2 )
eel.test2( y1, y2 )
}

\keyword{ Multivariate hypothesis testing }
\keyword{ non parametric test }
