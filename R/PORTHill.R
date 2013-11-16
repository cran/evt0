PORTHill = function (X,k =1,q =0.1) 
  
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
  
  osx = sort(X) 
        
  sum = 0 
  
  for (j in 1:k) 
  { 
    sum = sum + ( log(osx[n-j+1] - osx[nq] ) - log(osx[n-k] - osx[nq]) )  
     return(sum/k) 
  }
  
  
}