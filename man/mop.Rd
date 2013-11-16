\name{mop}
\alias{mop}

\title{
Mean of order p (MOP) statistics for the extreme value index (EVI) parameter
}
\description{
Mean of order p (MOP) statistics for the extreme value index (EVI) parameter.
}
\usage{
mop(X, k =1, p =0)
}

\arguments{
  \item{X}{ Data vector.
}
  \item{k}{ intermediate sequence, default is \code{k=1}.

}
  \item{p}{ mean order, default is \code{p=0}.
}
}
\details{
This method performs a mean of order p (is real) (MOP) methodology for estimating extreme value index (EVI) for heavy tailed models. At \code{p=0} MOP estimator is equal to classical Hill estimator.
}

\value{
 a numeric value
}
\references{
M.F. Brilhante, M.Ivette Gomes and Dinis Pestana. (2013). A simple generalisation of the Hill estimator
\emph{Computational Statistics and Data Analysis}, \bold{57}, 518-- 535.

J.Beran, D. Schell and M. Stehlik. (2013). The harmonic moment tail index estimator: asymptotic distribution and robustness  \emph{Ann Inst Stat Math},  Published Online.
}
\author{
B G Manjunath \email{bgmanjunath@gmail.com}
}

\examples{
set.seed(1.2345) 

# generate random samples               
X = rfrechet(500, loc = 0, scale = 1,shape = 1/0.5)

# estimate EVI 
mop(X,100, 0)
}

\keyword{mean of order p}
\keyword{Hill estimator}
