other.EVI = function (X,k =1, method=c("MO", "GH","MM")) 
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
  
 
  osx = sort(X) 
  
  method = match.arg(method)
  
  
  if (method == "MO") 
  {
    EVI.est = mo(osx,k)
  }
  
  else if (method == "GH")
  {
    EVI.est = gh(X,k)
  }  
  
  else if (method == "MM")
  {
    EVI.est = mm(osx,k)
  }  
  
  return(EVI.est)
  
  
}


# Moment(Mo) estimator
mo = function(osx,k)
{
  
  n = length(osx)
  M = length(2)
  
  for(i in 1:2)
  {
    sum = 0
    for (j in 1:k)
      sum = sum + ( log (osx[n-j+1]) - log (osx[n-k]) )^i
    M[i] = sum/k
    
  }
    
  moment = M[1] + (1/2)*( 1 -  ( (M[2]/M[1])^2 - 1)^(-1) )  
  
  return(moment)
  
}  


# Generalised Hill (GH) estimator
gh = function(X,k)
{
  
  evi =  mop(X,k,p=0,"MOP")
  
  sum =  0
  for (j in 1:k)
  {
    sum = sum + ( log(mop(X,j,p=0,"MOP")) - log(evi) )
  }
  
  return(evi + sum/k)
  
} 

# Mixed Moment(MM) estimator
mm = function(osx,k)
{
  
  n = length(osx)
  
  
  sum = 0
  for (j in 1:k)
      sum = sum + ( log (osx[n-j+1]) - log (osx[n-k]) )
  M = sum/k
  
  
  sum = 0
  for (j in 1:k)
    sum = sum + ( 1 - (osx[n-k]/osx[n-j+1])  )
  L = 1- sum/k
  
  phi = (M -L)/(L^2)
  
  return((phi-1)/(1+2*min(phi-1,0)))
  
}   

