\name{other.EVI}
\alias{other.EVI}
\alias{mo}
\alias{gh}
\alias{mm}
\title{
Other extreme value index estimates 
}
\description{
This function computes moment (MO), generalized Hill (GH) and mixed moment (MM) estimates for extreme value index (EVI).
}
\usage{
other.EVI(x, k, method = c("MO", "GH", "MM"))
}

\arguments{
  \item{x}{
Data vector.
}
  \item{k}{
a vector of number of upper order statistics.
}
  \item{method}{
Method used, moment estimate("MO", default), generalized Hill ("GH") and mixed moment ("MM").
}
}
\details{
Computation of moment and generalized Hill and mixed moment EVI estimators are based on the articles by Dekkers et al. (1989), Beirlant et al. (1996) and Fraga Alves et al. (2009), respectively.
}
\value{
a \code{k} dimensional vector of EVI estimates.
}
\references{
Dekkers, A., Einmahl, J. and L. de Haan. (1989). A moment estimator for the index of an extreme-value distribution. 
\emph{Ann. Statist.}, \bold{17}, 1833-- 1855.

Beirlant, J.,  Vynckier, P. and Teugels, J. (1996). Excess functions and estimation of the extreme-value index. 
\emph{Bernoulli}, \bold{2}, 293--318.

Fraga Alves, M.I.,  Gomes, M.I., de Haan, L. and Neves, C. (2009). The mixed moment estimator and location invariant alternatives. 
\emph{Extremes}, \bold{12}, 149--185.
}
\author{
B G Manjunath \email{bgmanjunath@gmail.com}, Frederico Caeiro \email{fac@fct.unl.pt}
}

\examples{
# generate random samples               
x = rfrechet(50000, loc = 0, scale = 1,shape = 1/0.5)

# estimate EVI 
other.EVI(x,c(500,5000,40000),"MO")

}

\keyword{moment}
\keyword{generalized Hill}
\keyword{mixed moment}
