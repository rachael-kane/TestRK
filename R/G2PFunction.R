G2P=function(X,h2,alpha,NQTN,distribution){
  n=nrow(X)
  m=ncol(X)
  #Sampling QTN
  QTN.position=sample(m,NQTN,replace=F)
  SNPQ=as.matrix(X[,QTN.position])
  QTN.position
  #QTN effects
  if(distribution=="rnorm")
  {addeffect=rnorm(NQTN,0,1)
  }else
  {addeffect=alpha^(1:NQTN)}
  #Simulate phenotype
  effect=SNPQ%*%addeffect
  effectvar=var(effect)
  residualvar=(effectvar-h2*effectvar)/h2
  residual=rnorm(n,0,sqrt(residualvar))
  y=effect+residual
  return(list(addeffect = addeffect, y=y, add = effect, residual = residual, QTN.position=QTN.position, SNPQ=SNPQ))
}
