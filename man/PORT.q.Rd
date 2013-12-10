\name{PORT.q}
\alias{PORT.q}
\title{
Peaks over random threshold high quantile estimate
}
\description{
This function computes high quantile or value-at-risk (VaR) estimate based on peaks over random threshold (PORT) Hill extreme value index (EVI) estimate.
}
\usage{
PORT.q(x, k, q1, q2, method = c("PMOP", "PRBMOP"))
}

\arguments{
  \item{x}{
Data vector.
}
  \item{k}{
a vector of number of upper order statistics.
}
  \item{q1}{
quantile for PORT.
}
  \item{q2}{
quantile level.
}
  \item{method}{
Method used, ("PMOP", default) 
                  and reduced-bias PMOP ("PRBMOP").
}
}
\details{
The computation of the high quantile estimate is based on the work by Gomes et al. (2006).
}
\value{
a \code{k} dimensional vector of PORT Hill and high quantile estimates. When \code{Method = "RBMOP"} shape and scale second order parameters estimates are also returned.
}
\references{
Araujo Santos, P., Fraga Alves, M.I. and Gomes, M.I. (2006). Peaks over random threshold methodology for tail index and quantile estimation.
\emph{Revstat}, \bold{4}(3), 227--247.

Gomes, M.I., Figueiredo, F., Henriques-Rodrigues, L. and Miranda, M.C. (2006). A quasi-PORT methodology for VaR based on second-order reduced-bias estimation.

Weissman, I. (1978). Estimation of parameters and large quantiles based on the k largest observations. \emph{J. Amer. Statist. Assoc.}, \bold{73}, 812-- 815.
}
\author{
B G Manjunath \email{bgmanjunath@gmail.com}
}

\seealso{
\code{\link[evt0]{PORT.Hill}}
}
\examples{
# generate random samples               
x = rfrechet(50000, loc = 0, scale = 1,shape = 1/0.5)

# estimate PORT Hill and quantile at level q2
PORT.q(x,c(1,500,5000),0.1,0.5,"PRBMOP")
}
\keyword{quasi-reduced bias}
\keyword{PORT}
\keyword{Hill}
\keyword{VaR}