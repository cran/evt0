\name{otherPORT.q}
\alias{otherPORT.q}
\title{
Other peaks over random threshold high quantile estimate
}
\description{
This function computes peaks over random threshold (PORT) high quantile or value-at-risk (VaR) based on  moment (MO), generalized Hill (GH) and mixed moment (MM) extreme value index (EVI) estimates. 
}
\usage{
otherPORT.q(x, k, q1, q2, method = c("MO", "GH", "MM"))
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
Method used, moment estimate("MO", default), generalized Hill ("GH") and mixed moment ("MM").
}
}
\details{
The computation of high quantile estimate is based on the method by Weissman (1978) and the EVI estimators are given in  Dekkers et al. (1989), Beirlant et al. (1996) and Fraga Alves et al. (2009).
}
\value{
a \code{k} dimensional vector of PORT EVI and high quantil estimates.
}
\references{
Araujo Santos, P., Fraga Alves, M.I. and Gomes, M.I. (2006). Peaks over random threshold methodology for tail index and quantile estimation.
\emph{Revstat}, \bold{4}(3), 227--247.

Dekkers, A., Einmahl, J. and L. de Haan. (1989). A moment estimator for the index of an extreme-value distribution. 
\emph{Ann. Statist.}, \bold{17}, 1833-- 1855.

Beirlant, J.,  Vynckier, P. and Teugels, J. (1996). Excess functions and estimation of the extreme-value index. 
\emph{Bernoulli}, \bold{2}, 293--318.

Fraga Alves, M.I.,  Gomes, M.I., de Haan, L. and Neves, C. (2009). The mixed moment estimator and location invariant alternatives. 
\emph{Extremes}, \bold{12}, 149--185.

Weissman, I. (1978). Estimation of parameters and large quantiles based on the k largest observations. \emph{J. Amer. Statist. Assoc.}, \bold{73}, 812-- 815.
}
\author{
B G Manjunath \email{bgmanjunath@gmail.com}
}

\seealso{
\code{\link[evt0]{other.EVI}}
}
\examples{
# generate random samples               
x = rfrechet(50000, loc = 0, scale = 1,shape = 1/0.5)

# estimate PORT EVI and high quantile at level q2
otherPORT.q(x,c(500,5000),0.1,0.5,"MO")
}

\keyword{PORT}
\keyword{moment}
\keyword{generalized Hill}
\keyword{mixed moment}
\keyword{VaR}
