\name{S_P500}
\alias{S_P500}
\docType{data}
\title{
S&P500
}
\description{
Log-returns of S&P500  Index from 05-01-1960 untill 16-10-1987.
}
\usage{data(S_P500)}
\format{
  A data frame with 6984 observations on the following variable.
  \describe{
    \item{\code{returns}}{a numeric vector}
  }
}
\details{
Log-returns of S&P500  Index from 05-01-1960 untill 16-10-1987.
}

\examples{
data(S_P500)
str(S_P500)
plot(S_P500$returns)
}
\keyword{datasets}
