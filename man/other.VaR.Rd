\name{other.VaR}
\alias{other.VaR}
\title{
Other method of value at risk (VaR) computation
}
\description{
This function computes VaR using moment, generalized Hill and mixed moment method.
}
\usage{
other.VaR(X, k = 1, q, method = c("MO", "GH", "MM"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
Data vector.
}
  \item{k}{
intermediate sequence, default is \code{k=1}.
}
  \item{q}{
quantile.
}
  \item{method}{
Method used, moment estimate("MO", default), generalized Hill ("GH") and mixed moment ("MM").
}
}
\details{
This function computes VaR using moment, generalized Hill and mixed moment method.
}
\value{
a numeric value
}
\references{
M. Ivette Gomes, M.F. Brilhante and Dinis Pestana. (2013). A Mean-of-order-p Class of Value-at-Risk Estimates \emph{Submitted article}.
}
\author{
B G Manjunath \email{bgmanjunath@gmail.com}
}

\examples{
set.seed(1.2345) 

# generate random samples               
X = rfrechet(5000, loc = 0, scale = 1,shape = 1/0.5)

# estimate VaR 
other.VaR(X,500,0.5,"MO")
}

\keyword{moment estimate}
\keyword{generalized Hill}
\keyword{mixed moment}
\keyword{value at risk}