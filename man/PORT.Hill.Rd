\name{PORT.Hill}
\alias{PORT.Hill}
\alias{PORT.MOP}
\alias{PORT.RBMOP}
\title{
Peaks over random threshold (PORT) Hill estimator
}
\description{
Peaks over random threshold (PORT) Hill estimator.
}
\usage{
PORT.Hill(X, k =1, q =0.1, method = c("MOP", "RBMOP"))
}

\arguments{
  \item{X}{
Data vector.
}
  \item{k}{
 intermediate sequence, default is \code{k=1}.
}
  \item{q}{
quantile, default is \code{q=0.1}.
}
  \item{method}{
Method used, ("MOP", default) 
                  and reduced-bias MOP ("RBMOP").
}
}
\details{
This method performs a peaks over random threshold (PORT) methodology for estimating extreme value index (EVI) for heavy tailed models.
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
PORT.Hill(X, 500, 0.2,"MOP")
}

\keyword{PORT}
\keyword{Hill estimator}