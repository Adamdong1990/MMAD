#' Plot the GaF object
#'
#' @param x The GaF object, see \code{\link{GaFrailtyMM}}.
#' @param XLAB x label, default is "Time".
#' @param YLAB y label, default is "Cumulative hazard".
#' @param TYPE type value, default is "s".
#' @param LTY lty value for line, default is 1.
#' @param LWD line width, default is 1.
#' @param COL color parameter, default is gray(0).
#' @param digits The digits after the decimal point, default = 4.
#' @param ... Additional arguments
#'
#' @return the dataframe of "Time" and accumulative hazard \eqn{\Lambda}.
#' @export
#'
#' @examples
#' library(survival)
#' result <- GaFrailtyMM(Surv(time, status) ~ age + sex + cluster(id), data=kidney)
#' \dontrun{ plot.GaF(result) }
#'
plot.GaF <-function(x, XLAB="Time", YLAB="Cumulative hazard", TYPE="s", LTY=1, LWD=1, COL=gray(0), digits=4, ...)
{
  LA <- sort(x$Lambda)
  time <- sort(x$time)

  tab <- data.frame(Time=round(time,digits), Lambda=round( LA,digits) )
  tab <- unique(tab)
  rownames(tab) <- 1:nrow(tab)
  print(tab)

  # Add a figure frame
  # The meaning of type="n" is not to add any elements to the graph,
  # but only to draw the coordinate axis
  XLIM <- range(time[time!=Inf])
  YLIM <- range(LA[LA!=Inf])
  plot(XLIM, YLIM, type="n", xlab=XLAB, ylab=YLAB, ...)

  do.call("lines", c(list(x=time,y=LA), type=TYPE, lty=LTY, lwd=LWD, col=COL,...))

}
