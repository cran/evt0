\name{mop.VaR}
\alias{mop.VaR}
\title{
Value at Risk (VaR) using mean of order p (MOP) 
}
\description{
This function computes VaR using MOP method.
}
\usage{
mop.VaR(X, k = 1, p = 0, q, method = c("MOP", "RBMOP"))
}

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
  \item{q}{quantile.
}
  \item{method}{
Method used, ("MOP", default) 
                  and reduced-bias MOP ("RBMOP").
}
}
\details{
This function computes VaR using MOP method.
}
\value{
a numeric value
}
\references{
M.F. Brilhante, M.Ivette Gomes and Dinis Pestana. (2013). A simple generalisation of the Hill estimator
\emph{Computational Statistics and Data Analysis}, \bold{57}, 518-- 535.

M. Ivette Gomes, M.F. Brilhante and Dinis Pestana. (2013). New reduced-bias estimators of a positive extreme value index \emph{Submitted article}.

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
mop.VaR(X,500,0,0.5,"MOP")
}

\keyword{mean of order p}
\keyword{VaR}
