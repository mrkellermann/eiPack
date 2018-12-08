\name{densityplot}
\alias{densityplot}
\alias{densityplot.lambdaMD}
\alias{densityplot.lambdaReg}
\alias{densityplot.lambdaRegBayes}

\title{Density plots for population level parameters}

\description{Generates a density plot for population level quantities of
  interest output by \code{\link{lambda.MD}}, \code{\link{lambda.reg}},
  and \code{\link{lambda.reg.bayes}}.  For the Bayesian methods,
  \code{densityplot} plots the kernel density for the draws.  For the
  frequentist \code{\link{lambda.reg}} method, \code{densityplot} plots
the canonical Normal density conditional on the mean and standard error
output by \code{\link{lambda.reg}}.}

\usage{
\method{densityplot}{lambdaMD}(x, by = "column", col, xlim, ylim,
             main = "", sub = NULL, xlab, ylab,
             lty = par("lty"), lwd = par("lwd"), ...)
\method{densityplot}{lambdaRegBayes}(x, by = "column", col, xlim, ylim,
             main = "", sub = NULL, xlab, ylab,
             lty = par("lty"), lwd = par("lwd"), ...)
\method{densityplot}{lambdaReg}(x, by = "column", col, xlim, ylim,
             main = "", sub = NULL, xlab, ylab,
             lty = par("lty"), lwd = par("lwd"), ...)
}

\arguments{
\item{x}{output from \code{\link{lambda.MD}}, \code{\link{lambda.reg}},
or \code{\link{lambda.reg.bayes}}.}
\item{by}{character string (defaulting to \code{"column"}) specifying
whether to panel the density plot by \code{"row"} or \code{"column"}
marginal.}
\item{col}{an optional vector of colors, with length corresponding to
the number of marginals selected in \code{by}.  Defaults to
\code{rainbow}.}
\item{xlim,ylim}{optional limits for the x-axis and y-axis, passed to
\code{plot}.}
\item{main,sub}{optional title and subtitle, passed to \code{plot}.}
\item{xlab,ylab}{optional labels for the x- and y-axes, passed to
\code{plot}.}
\item{lty,lwd}{optional arguments for line type and line width, passed
to \code{lines} and \code{plot}.  If either \code{lty} or \code{lwd} are
vectors, it must correspond to the number of row or column marginals selected.}
\item{...}{additional arguments passed to \code{par}.}
}
\value{
  A plot with density lines for the selected margin (row or column).  
}
\seealso{\code{plot}, \code{segments}, \code{par}}

\author{
  Olivia Lau <\email{olivia.lau@post.harvard.edu}>
}

\keyword{hplot}
