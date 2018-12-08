densityplot <- function(x, by = "column", col,
                         xlim = c(0,1), ylim,
                         main = "", sub = NULL, xlab,
                         ylab, lty = par("lty"), 
                         lwd = par("lwd"), ...) {
  UseMethod("densityplot", x)
}
