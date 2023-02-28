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
#' CoxProfile(x, I, t, be, la)
CoxProfile <- function(x, I, t, be, la, n, q)
{
  # compute LA
  LA=vector(length=n)
  for(i in 1:n){
    LA[i]=sum(la*(t<=t[i]))
  }
  # compute la
  AE=exp(x%*%be)
  de=abs(x)/matrix(apply(abs(x),1,sum),n,q)
  sum_1=vector(length=n)
  sum_2=sum_3=matrix(NA,n,q)
  for(i in 1:n){
    sum_1[i]=sum(AE*(t>=t[i]))
    for(p in 1:q){
      sum_2[i,p]=sum(AE*x[,p]*(t>=t[i]))/sum_1[i]
      sum_3[i,p]=sum(AE*x[,p]*x[,p]*(t>=t[i])/de[,p])/sum_1[i]
    }
  }
  newla=I/sum_1
  # compute beta
  newbe=vector(length=q)
  for(p in 1:q){
    f1=sum(I*x[,p]-I*sum_2[,p])
    f2=-sum(I*sum_3[,p])
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
  return(list(newlfun=newlfun, newbe=newbe, newla=newla, newLA=newLA, eps=eps))
}
