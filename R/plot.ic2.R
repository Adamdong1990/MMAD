plot.ic2 <-function(object, XLAB="Time", YLAB="Survival", dtype="l", main="Survival Function",
                    LTY=1:9, LWD=1, col=gray(0),XLEG=NULL,YLEG=NULL, ...)
{
  #XLEG<-max(0,min(object$intmap[1,]))
  #YLEG<-.1

  XLIM <- object$s
  YLIM <- object$S

  plot(XLIM, YLIM, xlab=XLAB, ylab=YLAB, type=dtype, main=main, lty=LTY, lwd=LWD, col=col,...)
}
