rm(list=ls())

sam=function(n,be,th,la){
  q=length(be)
  x=matrix(rnorm(n*q,-1,1),n,q)
  u=runif(n)
  t=-log(u)/(la*exp(x%*%be))
  #t=(exp(-log(u)/exp(x%*%be))-1)/la
  cen=3.7
  I=1*(t<=cen)
  t=pmin(t,cen)
  return(list(x=x,I=I,t=t))
}

n=200;q=3
be.true=matrix(c(-0.5,1,2),ncol=1)
t.true=c(0.1,0.5)
LA.TRUE=2*t.true

be.in=be.true*0.5

da=sam(n,be=be.true,la=2)
data=da$x;time=da$t;status=da$I;be=be.in;bootstrap = TRUE

result = coxmm(data, time, status, be, convergence=1e-6, b=1000, method = "Pro")
summary.cox(result)

result = coxmm(data, time, status, be, convergence=1e-8, b=100, method = "NPro")
summary.cox(result)


t1=time[1:9];t2=time[10:n]
non.LA1=result$LA[1:9]
non.LA2=result$LA[10:n]

t1=sort(time)
non.LA1=sort(result$LA)

plot(t1,non.LA1,pch=22,col="red",xlab="month",ylab="estimated baseline cumulative hazard rate",cex.lab=0.8,cex.axis=0.8)

plot(t1,non.LA1,col="red",xlab="month",ylab="estimated baseline cumulative hazard rate")

points(t1,non.LA1,lty=6,type="s",col="black")

