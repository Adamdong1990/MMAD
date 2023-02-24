Aintmap<-function(L,R,Lin=NULL,Rin=NULL){
  n<-length(L)
  if (is.null(Lin) & is.null(Rin)){
    Lin<-rep(FALSE,n)
    Rin<-rep(TRUE,n)
    Lin[L==R]<-TRUE
    Rin[R==Inf]<-FALSE
  } else if (length(Lin)==1 & length(Rin)==1 & is.logical(Lin) & is.logical(Rin) ){
    Lin<-rep(Lin,n)
    Rin<-rep(Rin,n)
  } else if (length(Lin)!=n | length(Rin)!=n | !all(is.logical(Lin)) | !all(is.logical(Rin)) ){
    stop("Lin and Rin should be either NULL, logical length 1 or length same as L,R")
  }
  if(n != length(R))
    stop("length of L and R must be the same")
  # calculate a small number, eps, to differentiate between e.g.,  [L,R] and (L,R]
  # we will treat, (L,R] as [L+eps,R], and [L,R) as [L,R-eps]
  # since eps is only
  # used in ranking we do not need to make it super small
  # just smaller than the smallest difference
  LRvalues<-sort(unique(c(0,L,R,Inf)))
  eps<- min(diff(LRvalues))/2
  Le<-L
  Re<-R
  Le[!Lin]<-L[!Lin]+eps
  Re[!Rin]<-R[!Rin]-eps
  # let s be the vector of ordered L and R values with
  # R values later when there are ties
  # then intmap are values s[i] and s[i+1] where s[i] is
  # associated with L and s[i+1] is associated with R
  oLR<-order(c(Le,Re+eps/2) )
  # find the Turnbull intervals, or innermost intervals
  # this is the same as the primary reduction of
  ### Aragon and Eberly (1992) J of Computational and Graphical
  ###     Statistics 1:129-140
  # label L=1 and R=2
  Leq1.Req2<-c(rep(1,n),rep(2,n))
  # order and see if an R is followed by an L
  # take difference of Leq1.Req2 after putting them in
  # order, then if the difference is 1 then the R=2 is followed by L=1
  flag<- c(0,diff( Leq1.Req2[oLR] ))
  R.right.of.L<- (1:(2*n))[flag==1]
  intmapR<- c(L,R)[oLR][R.right.of.L]
  intmapL<- c(L,R)[oLR][R.right.of.L - 1]
  intmapRin<- c(Lin,Rin)[oLR][R.right.of.L]
  intmapLin<- c(Lin,Rin)[oLR][R.right.of.L - 1]
  intmap<-matrix(c(intmapL,intmapR),byrow=TRUE,nrow=2)
  attr(intmap,"LRin")<-matrix(c(intmapLin,intmapRin),byrow=TRUE,nrow=2)
  k<-dim(intmap)[[2]]
  Lbracket<-rep("(",k)
  Lbracket[intmapLin]<-"["
  Rbracket<-rep(")",k)
  Rbracket[intmapRin]<-"]"
  intname<-paste(Lbracket,intmapL,",",intmapR,Rbracket,sep="")
  A<-matrix(0,n,k,dimnames=list(1:n,intname))
  intmapLe<-intmapL
  intmapLe[!intmapLin]<-intmapL[!intmapLin]+eps
  intmapRe<-intmapR
  intmapRe[!intmapRin]<-intmapR[!intmapRin]-eps
  for (i in 1:n){
    tempint<- Le[i]<=intmapRe & Re[i]>intmapLe
    A[i,tempint]<-1
  }

  # previous versions (<=0.9-9.1) did primary reduction twice,
  # once as described in Turnbull (see above) and once as
  # described in Aragon and Eberly (1992, J of Computational and Graphical
  #    Statistics 1:129-140)
  # both do same thing, so we do not need to do it twice

  ## fix error when intmap=(0,Inf) and k=1, previously A was column matrix of 0, should be a column matrix of 1
  if (k==1 & intmap[1,1]==0 & intmap[2,1]==Inf) A[A==0]<-1

  out<-list(A=A,intmap=intmap,intmapL=intmapL,intmapR=intmapR)
  out
}
