\name{mop.AREFF}
\alias{mop.AREFF}
\title{
Asymptotic efficiency of mean of order p 
}
\description{
This function compute asymptotic relative efficiency of mean of order p (MOP) with respect to classical Hill estimator.
}
\usage{
mop.AREFF(x, k, p)
}

\arguments{
  \item{x}{
Data vector.
}
  \item{k}{
a vector of number of upper order statistics.
}
  \item{p}{
a vector of mean order.
}
  
}

\details{
Given two biased estimators MOP and Hill, the asymptotic root efficiency (AREFF) of MOP relatively to Hill is given in Brilhante et al. (2013). Note that highest the AREFF indicator the better is the MOP estimator.
}
\value{
a matrix of asymptotic relative efficiency estimates, corresponds to \code{k} row and \code{p} columns.
}

\references{
Brilhante, M.F., Gomes, M.I. and Pestana, D. (2013). A simple generalisation of the Hill estimator.
\emph{Computational Statistics and Data Analysis}, \bold{57}, 518-- 535.
}
\author{
B G Manjunath \email{bgmanjunath@gmail.com}
}

\seealso{
\code{\link[evt0]{mop}}
}
\examples{
# generate random samples               
x = rfrechet(50000, loc = 0, scale = 1,shape = 1)

# estimate AREFF
mop.AREFF(x,c(1,500,5000,49999), c(-1,0,0.1))
}

\keyword{MOP}
\keyword{AREFF}
