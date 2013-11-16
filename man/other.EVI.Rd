\name{other.EVI}
\alias{other.EVI}
\alias{mo}
\alias{gh}
\alias{mm}
\title{
Other extreme value index(EVI) estimate mathods
}
\description{
This function computes moment, generalized Hill and mixed moment estimates for EVI.
}
\usage{
other.EVI(X, k = 1, method = c("MO", "GH", "MM"))
}

\arguments{
  \item{X}{
Data vector.
}
  \item{k}{
intermediate sequence, default is \code{k=1}.
}
  \item{method}{
Method used, moment estimate("MO", default), generalized Hill ("GH") and mixed moment ("MM").
}
}
\details{
This function computes moment, generalized Hill and mixed moment estimates for EVI.
}
\value{
a numeric value
}
\references{
M.F. Brilhante, M.Ivette Gomes and Dinis Pestana. (2013). A simple generalisation of the Hill estimator
\emph{Computational Statistics and Data Analysis}, \bold{57}, 518-- 535.

}
\author{
B G Manjunath \email{bgmanjunath@gmail.com}
}

\examples{
set.seed(1.2345) 

# generate random samples               
X = rfrechet(5000, loc = 0, scale = 1,shape = 1/0.5)

# estimate VaR 
other.EVI(X,500,"MO")
}

\keyword{moment estimate}
\keyword{generalized Hill}
\keyword{mixed moment}

