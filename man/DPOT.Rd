\name{DPOT}
\alias{DPOT}
\title{
Duration based peaks over threshold  value-at-risk  forecast
}
\description{
This function calculate the value-at-risk (VaR) forecast for the durations-based peaks over threshold (DPOT) models.
}
\usage{
DPOT(x, cov=0.01, c=0.75, th=0.1,  nd=1000)
}

\arguments{
  \item{x}{
Data vector.
}
  \item{cov}{
Coverage value, default is \code{cov=0.01}.
}
  \item{c}{
Tuning parameter, default is \code{c=0.75}.
}
  \item{th}{
Threshold value, default is \code{th=0.1}.
}
  \item{nd}{
Returns days, default is \code{nd=1000}.
}
}
\details{
In financial time series a relation between the excesses and the durations between excesses is usuallly observed.  Araujo Santos and Fraga Alves (2013) propose using this dependece to improve the risk forecasts with DPOT models. The computation method in \code{DPOT()} function is based on the work from Araujo Santos and Fraga Alves (2012).
}

\section{Warning}{
After running the function following message appears:
\code{In log(1+gamma*y/(alpha1*(1/x)^c )): NaNs produced} when the gamma is negative but the optimizer continue to other iternations choosing other values until it converge.
}

\value{
VaR forecast and also MLE estimates of shape and time scale parameters.
}

\references{

Araujo Santos, P. and Fraga Alves, M.I. (2013). Forecasting Value-at-Risk with a duration-based POT method. \emph{Mathematics and Computers in Simulation}, \bold{94}, 295--309.

Araujo Santos, P. and Fraga Alves, M.I. (2012). R Program to Implement the DPOT Model. \emph{Unpublished article}.

}
\author{
P. Araujo Santos \email{paulo.santos@esg.ipsantarem.pt}, M.I. Fraga Alves \email{isabel.alves@fc.ul.pt}
}

\examples{
#Read S&P500 from data file
data(S_P500)
str(S_P500)

# One day ahead VaR forecast
DPOT(S_P500$returns,0.01,0.75,0.1,1000)
}

\keyword{ DPOT}
\keyword{ VaR }
