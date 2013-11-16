\name{PORT.VaR}
\alias{PORT.VaR}
\title{
Peaks over random threshold (PORT) value at risk  (VaR)
}
\description{
This function computes VaR using PORT method.
}
\usage{
PORT.VaR(X, k = 1, q1 =0.1, q2, method = c("MOP", "RBMOP"))
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
Method used, ("MOP", default) 
                  and reduced-bias MOP ("RBMOP").
}
}
\details{
This function computes VaR using PORT method.
}
\value{
a numeric value
}
\references{
P.A. Santos, M. Isabel Fraga Alves and M. Ivette Gomes. (2006). Peaks over random threshold methodology for tail index and quantile estimation.
\emph{Revstat}, \bold{4}(3), 227--247.
}
\author{
B G Manjunath  \email{bgmanjunath@gmail.com}
}

\examples{
set.seed(1.2345) 
              
# generate random samples               
X = rfrechet(5000, loc = 0, scale = 1,shape = 1/0.5)

# estimate EVI 
PORT.VaR(X, 500, 0.2,0.5,"MOP")
}

\keyword{PORT}
\keyword{VaR}