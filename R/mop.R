mop = function (X,k =1,p =0, method=c("MOP", "RBMOP")) 
  
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
    stop("Some of the data points are not real.")
  }
  
  if (!is.numeric(p)) 
  {
    stop("p must be a real.")
  }
  
  osx = sort(X) 
  
  method = match.arg(method)
  
  
  if (method == "MOP") 
  {
    EVI.est = mop.MOP(osx,k,p)
  }
  
  else if (method == "RBMOP")
  {
    EVI.est = mop.RBMOP(osx,k,p)
  }  
  
  
  return(EVI.est)  
  
}

# mean of order p extreme value index estimate

mop.MOP = function(osx,k,p)
{
  
  n = length(osx)
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

# Reduced bias mean of order p extreme value index estimate

mop.RBMOP = function(osx,k,p)
{ 
  n = length(osx)
  
  H =  mop.MOP(osx,k,p)
  
  rhoest  = mop.rho(osx)
  
  betaest = mop.beta(osx,rhoest)
  
  estimate = H* (1 -   (  (betaest* (1-p*H)) / (1-rhoest-p*H)  )* (n/k)^(rhoest)    )
  
  return(list(EVI =estimate, rho =rhoest,beta=betaest))
}  

# Rho estimation
mop.rho = function(osx)
{
  
  n = length(osx)
  M = numeric(3)
  krho = floor(n^(0.995)):floor(n^(0.999))
  
  tau0 = numeric(length(krho))
  tau1 = numeric(length(krho))                              
  
  for (l in 1:length(krho))
  {
    K = krho[l]
    
    for (i in 1:3)
    {
      sum =  0
      for (j in 1:K)
        sum = sum + (  (osx[n-j+1]) / (osx[n-K]) )^i     
      M[i] =  sum/K        
    }
    
    W0 =  ( log(M[1]) - (1/2)* log(M[2]/2) )/ ((1/2)* log(M[2]/2) - (1/3)* log(M[3]/6))
    
    W1 =  ( M[1] - (M[2]/2)^(1/2)  )/ ( (M[2]/2)^(1/2) -  (M[3]/6)^(1/3) )
    
    tau0[l] = -abs( 3*(W0-1)/(W0 -3) )
    tau1[l]= -abs( 3*(W1-1)/(W1 -3) )    
  }
  
  
  
  tau0median = median(tau0)
  tau1median = median(tau1)
  
  
  
  sumtau0 = as.numeric((tau0 - tau0median) %*% (tau0 - tau0median))
  sumtau1 = as.numeric((tau1 - tau0median) %*% (tau1 - tau0median))
  
  if((sumtau0 < sumtau1) || (sumtau0 == sumtau1) )
  {
    return(tau0[l])   
    
  }
  
  else
  {
    return(tau1[l])
  }
    
  
}

# Beta estimation

mop.beta = function(osx,rhoest)
{
  n = length(osx)
  
  k1 = floor(n^(0.999))
  
  v = c(0,rhoest,2*rhoest)
  D = numeric(length(v))
  
  
  for (i in 1:length(v))
  {
    sum =  0
    for (j in 1:k1)
      sum = sum + ( (j/k1) )^(-v[i]) * j*(  (osx[n-j+1]) / (osx[n-k1]) )
    D[i]= sum/k1
  }
  
  
  sum =  0
  for (j in 1:k1)
    sum = sum + ( (j/k1) )^(-v[2]) 
  d= sum/k1
  
  
  const =  (d*D[0] -D[1])/(d*D[1]-D[2]) 
  
  betaest  =  (k1/n)^(v[2])*(d*D[1] -D[2])/(d*D[2]-D[3]) 
  return(betaest)
  
  
}

