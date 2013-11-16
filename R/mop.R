 mop = function (X,k =1,p =0) 
  
{ 
    
  n = length(X) 
  
  if (n  < 2) {
    stop("Data vector X with at least two sample points required.")
  }
  
    
  if (k < 1 || k > n || k != as.integer(k))
  {
    stop("k must be greater than or equaul to 1 and less than sample size.")
  }
    
  if (!is.numeric(X)) {
    stop("some of the data points are not real.")
  }
  
  if (!is.numeric(p)) 
  {
    stop("p must be a  real.")
  }
  
  
     
  osx = sort(X) 
  
  sum = 0
  if (p == 0)
  {
      
      for (j in 1:k)
          sum = sum + ( log (osx[n-j+1]) - log (osx[n-k]) )
          return(sum / k)
      
  }
    
  else 
  {
        
      for (j in 1:k)
        sum = sum + (  (osx[n-j+1]) / (osx[n-k]) )^p      
        return((1 - (sum/k)^-1) / p) 
      
  } 
    
 
  
}