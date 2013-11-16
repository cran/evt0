mop.AREFF = function (X,k =1,p =0) 
  
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
  
  if (!is.numeric(p)) 
  {
    stop("p must be a  real.")
  }
  
   
  
  EVI.est = mop(X,k,p,"MOP")
  
  RBEVI   = mop(X,k,p,"RBMOP")
  
  rhoest = as.numeric(RBEVI[2])

 if(p*EVI.est>0.5)
 {
    stop("Condition on EVI and p  is not satified.")
  }
  
  c1  =  ( sqrt(1-2*p*EVI.est)/(1-p*EVI.est)) ^(-2*rhoest)

  c2 = abs( (1-p*EVI.est-rhoest)/((1-rhoest)*(1-p*rhoest)) )
  
  return( (c1*c2)^(1/(1-2*rhoest)) )
  
  
}