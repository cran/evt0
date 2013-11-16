other.VaR = function (X,k=1,q, method=c("MO", "GH","MM"))   
{ 
  
  n = length(X) 
  
  if (n  < 2) {
    stop("Data vector X with at least two sample points required.")
  }
  
  
  if (k < 1 || k > n || k != as.integer(k))
  {
    stop("k must be greater than or equal to 1 and less than sample size.")
  }
  
  if (!is.numeric(X)) {
    stop("some of the data points are not real.")
  }
  
    
  if ( q < 0 || q > 1 || !is.numeric(q)) 
  {
    stop("q must lie between 0 and 1.")
  }
  
  method = match.arg(method)
  
  
  EVI.est = other.EVI(X,k,method)
  
  osx = sort(X) 
  
  
  otherVaR = (osx[n-k]) *  (k/(n*q))^EVI.est
  
  return(otherVaR) 
  
}
