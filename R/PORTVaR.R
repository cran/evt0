PORT.VaR = function (X,k=1,q1=0.1,q2,method=c("MOP", "RBMOP")) 
  
{ 
   n = length(X) 
  
  nq = floor(n*q1) + 1
  
  if (nq < 2 ) 
  {
    stop("Insufficient data vector, check sample size and q.")
  }
  
  
  if (k < 1 || k > n-nq -1  || !is.numeric(k)   || k != as.integer(k)  )
  {
    stop("k must be greater than or equal to 1 and less than exceedance sample size.")
  }
  
  if (!is.numeric(X)) {
    stop("some of the data points are not real.")
  }
  
  
  if ( q1 < 0 || q1 > 1 || !is.numeric(q1)) 
  {
    stop("q1 must lie between 0 and 1.")
  }
  
  
  if ( q2 < 0 || q2 > 1 || !is.numeric(q2)) 
  {
    stop("q2 must lie between 0 and 1.")
  }
  
   
  method = match.arg(method)
  
  osx = sort(X) 
   
    
  PORT.EVI = PORT.Hill(X,k,q1,method)
   
  
  if (method == "RBMOP") 
  {
    PORT.EVI = as.numeric(PORT.EVI[1])
  }
  
  PORTVaR =   (osx[n-k]- osx[nq] )  * (k/(n*q2))^PORT.EVI + osx[nq]
  
  return(PORTVaR)
  
  
}