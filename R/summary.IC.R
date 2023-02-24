summary.IC<-function(object,digits=4,epsilon=4,...)
{
  error <- object$print_err
  convergence <- object$convergence
  if (error > convergence)
  {
    cat("Convergence result: Did Not Converge")
    cat("\n")
  }
  else
  {
    p <- object$p
    intmap<-object$intmap
    k<-dim(intmap)[[2]]

    # make interval description match LRin attributes
    # ( and ) denote LRin value is FALSE
    # [ and ] denote LRin value is TRUE
    LRin<- attr(intmap,"LRin")
    # default to [ and ] if LRin is not given
    # if (is.null(LRin)) LRin<-matrix(TRUE,2,k)
    Lbracket<-rep("(",k)
    Lbracket[LRin[1,]]<-"["
    Rbracket<-rep(")",k)
    Rbracket[LRin[2,]]<-"]"
    intname<-paste(Lbracket,round(intmap[1,],digits),",",round(intmap[2,],digits),Rbracket,sep="")

    tab <- data.frame(Interval=intname, Probability=round(p,epsilon))
    df_tab <- tab[which(tab$Probability>0),]
    # rename
    row.names(df_tab) <- c(1:length(df_tab$Probability))

    print(df_tab)
  }
}
