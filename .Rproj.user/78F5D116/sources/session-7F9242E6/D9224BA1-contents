#' Title
#'
#' @param x a data frame containing the variables of the model
#' @param I The status indicator, normally 0=alive, 1=dead.
#' @param t for right censored data
#' @param be a vector of regression parameters
#' @param la baseline hazard rate
#'
#' @return requires a value
#' @export
#'
#' @examples
#' x = matrix(rnorm(10*1,-1,1),10,1)
#' la = 2; be = 0.5
#' u=runif(10); t=-log(u)/(la*exp(x%*%be))
#' cen=3.7; I=1*(t<=cen)
#' Profile(x, I, t, be, la)
Profile <- function(x, I, t, be, la)
{
  n=dim(x)[1]
  q=dim(x)[2]
  # compute LA
  LA=vector(length=n)
  for(i in 1:n){
    LA[i]=sum(la*(t<=t[i]))
  }
  # compute lambda
  AE=exp(x%*%be)
  newla=vector(length=n)
  for(i in 1:n){
    newla[i]=I[i]/sum(AE*(t>=t[i]))
  }
  # compute beta
  b=LA*exp(x%*%be)
  de=abs(x)/matrix(apply(abs(x),1,sum),n,q)
  newbe=vector(length=q)
  for(p in 1:q){
    f1=sum(I*x[,p]-x[,p]*b)
    f2=-sum(2*x[,p]*x[,p]*b/de[,p])
    newbe[p]=be[p]-f1/f2
  }
  # first
  a1=(la*exp(x%*%be))^I
  a2=LA*exp(x%*%be)
  lfun=sum(log(a1)-a2)
  # second
  newLA=vector(length=n)
  for(i in 1:n){
    newLA[i]=sum(newla*(t<=t[i]))
  }
  a1=(newla*exp(x%*%newbe))^I
  a2=newLA*exp(x%*%newbe)
  newlfun=sum(log(a1)-a2)
  #########################
  eps=abs(newlfun-lfun)/(abs(lfun)+1)
  return(list(newbe=newbe, newla=newla, newLA=newLA, eps=eps))
}
