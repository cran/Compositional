\name{The Helmert sub-matrix}
\alias{helm}
\title{
The Helmert sub-matrix
}

\description{
The Helmert sub-matrix.
}

\usage{
helm(n)
}

\arguments{
  \item{n}{
  A number grater than or equal to 2.
}
}

\details{
The Helmert sub-matrix is returned. It is an orthogonal matrix without the first row.
}
\value{
A \eqn{(n-1) \times n} matrix.
}

\references{
Tsagris M.T., Preston S. and Wood A.T.A. (2011). A data-based power transformation for compositional data.
In Proceedings of the 4th Compositional Data Analysis Workshop, Girona, Spain. http://arxiv.org/pdf/1106.1451.pdf

John Aitchison (2003). The Statistical Analysis of Compositional Data p. 99 Blackburn Press.

Lancaster H. O. (1965). The Helmert matrices. The American Mathematical Monthly 72(1): 4-12.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <gioathineou@gmail.com>.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{alfa}, \link{alfainv}
}
}

\examples{
helm(3)
helm(5)
}

\keyword{ orthogonal matrix }
