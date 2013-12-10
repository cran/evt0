\name{mop.q}
\alias{mop.q}
\title{
High qunatile estimate  by mean of order p  statistic
}
\description{
This function compute estimate of high quantile or value-at-risk (VAR) using mean of order p (MOP) method.
}
\usage{
mop.q(x, k, p, q, method = c("MOP", "RBMOP"))
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
  \item{q}{
quantile level.
}
  \item{method}{
  Method used, ("MOP", default) 
                  and reduced-bias MOP ("RBMOP").

}
}
\details{
For heavy tails, Gomes et al. (2013) introduces a new class of high quantile estimators based on a class of mean of order p (MOP) extreme value index (EVI) estimators is givin by
\deqn{Q(k) = (X_{n-k:n}) (k/nq)^{H_p(k)},}
where \eqn{H_p(k)} is MOP EVI estimator and \eqn{X_{i:n}} is order statistic.
}
\value{
a matrix of EVI and VaR estimates, corresponds to \code{k} row and \code{p} columns. When \code{Method = "RBMOP"} shape and scale second order parameters estimates are also returned.
}
\references{
Brilhante, M.F., Gomes, M.I. and Pestana, D. (2013). A simple generalisation of the Hill estimator.
\emph{Computational Statistics and Data Analysis}, \bold{57}, 518-- 535.

Beran, J., Schell, D. and Stehlik, M. (2013). The harmonic moment tail index estimator: asymptotic distribution and robustness.  \emph{Ann Inst Stat Math},  Published Online.

Gomes, M.I., Brilhante, M.F. and Pestana, D. (2013). New reduced-bias estimators of a positive extreme value index. \emph{Submitted article}.

Weissman, I. (1978). Estimation of parameters and large quantiles based on the k largest observations. \emph{J. Amer. Statist. Assoc.}, \bold{73}, 812-- 815.
}

\author{
B G Manjunath \email{bgmanjunath@gmail.com}
}

\seealso{
\code{\link[evt0]{mop}}
}
\examples{
# generate random samples               
x = rfrechet(50000, loc = 0, scale = 1,shape = 1/0.5)

# estimate EVI and high quantile at level q
mop.q(x,c(1,500,5000,49999), c(-1,0,1),0.5,"RBMOP")

}

\keyword{ MOP }
\keyword{ EVI }
\keyword{reduced-bias}
\keyword{Hill}
\keyword{VaR}
