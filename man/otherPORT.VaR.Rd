\name{otherPORT.VaR}
\alias{otherPORT.VaR}
\title{
Other peaks over random threshold (PORT) value at risk (VaR) computation
}
\description{
This function compute PORT VaR using  moment, generalized Hill and mixed moment methodolgy. 
}
\usage{
otherPORT.VaR(X, k, q1=0.1, q2, method = c("MO", "GH", "MM"))
}

\arguments{
  \item{X}{
Data vector.
}
  \item{k}{
intermediate sequence, default is \code{k=1}.
}
  \item{q1}{
quantile, default is \code{q=0.1}.
}
  \item{q2}{
quantile for VaR
}
  \item{method}{
Method used, moment estimate("MO", default), generalized Hill ("GH") and mixed moment ("MM").
}
}
\details{
This function compute PORT VaR using  moment, generalized Hill and mixed moment methodolgy. 
}
\value{
a numeric value
}
\references{
M. Ivette Gomes, M.F. Brilhante and Dinis Pestana. (2013). A Mean-of-order-p Class of Value-at-Risk Estimates \emph{Submitted article}.
}
\author{
B G Manjunath  \email{bgmanjunath@gmail.com}
}

\examples{
set.seed(1.2345) 
              
# generate random samples               
X = rfrechet(5000, loc = 0, scale = 1,shape = 1/0.5)

# estimate EVI 
otherPORT.VaR(X, 500, 0.2,0.5,"MO")
}

\keyword{PORT}
\keyword{VaR}