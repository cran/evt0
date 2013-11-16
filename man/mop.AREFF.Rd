\name{mop.AREFF}
\alias{mop.AREFF}
\title{
Asymptotic efficiency of mean of order p (MOP) 
}
\description{
This function compute asymptotic relative efficiency of MOP with respect to classical Hill estimator.
}
\usage{
mop.AREFF(X, k = 1, p = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
Data vector.
}
  \item{k}{
  intermediate sequence, default is \code{k=1}.
}
  \item{p}{
mean order, default is \code{p=0}.
}
}
\details{
This function compute asymptotic relative efficiency of MOP with respect to classical Hill estimator.
}
\value{
a numeric value
}
\references{
M.F. Brilhante, M.Ivette Gomes and Dinis Pestana. (2013). A simple generalisation of the Hill estimator
\emph{Computational Statistics and Data Analysis}, \bold{57}, 518-- 535.
}
\author{
B G Manjunath \email{bgmanjunath@gmail.com}
}

\seealso{
\code{\link{mop}}
}
\examples{
set.seed(1.2345) 

# generate random samples               
X = rfrechet(5000, loc = 0, scale = 1,shape = 1/0.5)

# Asymptotic relative efficiency
mop.AREFF(X,100, 0)
}

\keyword{ mean of order P }
\keyword{ asymptotic efficiency}% __ONLY ONE__ keyword per line
