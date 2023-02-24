summary.IC2<-function(object,digits=4,epsilon=4,...)
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
    s <- round(object$s,digits)
    k <- length(s)-1
    Lbracket <- rep("(",k)
    Rbracket <- rep("]",k)
    intmapL <- s[1:k]
    intmapR <- s[1:k+1]
    intname<-paste(Lbracket,intmapL,",",intmapR,Rbracket,sep="")

    tab <- data.frame(Interval=intname, Probability=round(p[1:k+1],epsilon))
    df_tab <- tab[which(tab$Probability>0),]
    # rename
    row.names(df_tab) <- c(1:length(df_tab$Probability))

    print(df_tab)
  }
}
