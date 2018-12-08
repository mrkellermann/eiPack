\name{plot.bounds}
\alias{plot}
\alias{plot.bounds}
\title{Plot of deterministic bounds for units satisfying row thresholds}

\description{Plots the deterministic bounds on the proportion of
  row members within a specified column.}

\usage{
\method{plot}{bounds}(x, row, column, labels = TRUE, order = NULL,
    intersection = TRUE, xlab, ylab, col = par("fg"), 
    lty = par("lty"), lwd = par("lwd"), ...)
}

\arguments{
  \item{x}{output from \code{\link{bounds}}}
  \item{row}{a character string specifying the row of interest}
  \item{column}{a character string specifying the column of interest}
  \item{labels}{a logical toggle specifying whether precinct labels should
    be printed above interval bounds}
  \item{order}{an optional vector of values between 0 and 1 specifying
    the order (left-to-right) in which interval bounds are plotted}
  \item{intersection}{a logical toggle specifying whether the
    intersection of all plotted bounds (if it exists) should be plotted}
  \item{xlab, ylab, ...}{additional arguments passed to \code{plot}}
  \item{col, lty, lwd}{additional arguments passed to \code{segments}}
}
\value{
  A plot with vertical intervals indicating the deterministic bounds on
  the quantity of interest, and (optionally) a single horizontal
  interval indicating the intersection of these unit bounds.
}
\seealso{\code{bounds}}

\author{
  Ryan T. Moore <\email{rtm@american.edu}>
}

\keyword{hplot}
