\name{PORT.Hill}
\alias{PORT.Hill}
\alias{PORT.MOP}
\alias{PORT.RBMOP}
\title{
Peaks over random threshold  Hill estimate
}
\description{
This function performs peaks over random threshold (PORT) Hill methodology for estimating extreme value index (EVI) for heavy tailed models.
}
\usage{
PORT.Hill(x, k, q, method = c("PMOP", "PRBMOP"))
}

\arguments{
  \item{x}{
Data vector.
}
  \item{k}{
a vector of number of upper order statistics.
}
  \item{q}{
quantile for PORT.
}
  \item{method}{
Method used, ("PMOP", default) 
                  and reduced-bias PMOP ("PRBMOP").
}
}
\details{
The computation of PORT Hill estimator is based on the work by Araujo Santos et al. (2006). Reduced biased PORT Hill computation is based on quasi-PORT methodology, see Gomes et al.
}

\value{
a \code{k} dimensional vector of PORT Hill estimates. When \code{Method = "RBMOP"} shape and scale second order parameters estimates are also returned.
}

\references{
Araujo Santos, P., Fraga Alves, M.I. and Gomes, M.I. (2006). Peaks over random threshold methodology for tail index and quantile estimation.
\emph{Revstat}, \bold{4}(3), 227--247.

Gomes, M.I., Figueiredo, F., Henriques-Rodrigues, L. and Miranda, M.C. (2006). A quasi-PORT methodology for VaR based on second-order reduced-bias estimation.
}
\author{
B G Manjunath \email{bgmanjunath@gmail.com}, Frederico Caeiro \email{fac@fct.unl.pt}
}


\examples{
# generate random samples               
x = rfrechet(50000, loc = 0, scale = 1,shape = 1/0.5)

# estimate PORT Hill 
PORT.Hill(x,c(1,500,5000),0.1,"PRBMOP")
}

\keyword{quasi-reduced bias}
\keyword{PORT}
\keyword{Hill}
