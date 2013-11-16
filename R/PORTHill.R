PORT.Hill = function (X,k=1,q=0.1,method=c("MOP", "RBMOP")) 
  
{ 
    
  n = length(X) 
  
  nq = floor(n*q) + 1
  
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
  
  
  if ( q < 0 || q > 1 || !is.numeric(q)) 
  {
    stop("q must lie between 0 and 1.")
  }
  
  
  method = match.arg(method)
  
  
  if (method == "MOP") 
  {
    PORT.EVI = PORT.MOP(X,k,q)
  }
  
  else if (method == "RBMOP")
  {
    PORT.EVI = PORT.RBMOP(X,k,q)
  }  
  
  
  return(PORT.EVI)  
  
   
}

# PORT Hill estimator

PORT.MOP = function(X,k,q)
{
  
  n = length(X)
  osx = sort(X) 
  nq = floor(n*q) + 1
  
  sum = 0 
  
  for (j in 1:k) 
  { 
    sum = sum + ( log(osx[n-j+1] - osx[nq] ) - log(osx[n-k] - osx[nq]) )  
    return(sum/k) 
  }
    
}


# quasi-PORT Hill estimator

PORT.RBMOP = function(X,k,q)
{ 
  n = length(X)
  osx = sort(X) 
  
  nq = floor(n*q) + 1
  
  H =  PORT.MOP(X,k,q)
  
  S =  mop(X,k,p=0,"RBMOP")
  
  Sest  = as.integer(S)
  
  
  estimate = H* (1 -   ( Sest[3] / (1-Sest[2])  )* (n/k)^(Sest[2])    )
  
  return(list(PORTEVI =estimate, rho =Sest[2],beta =Sest[3]))
}  
