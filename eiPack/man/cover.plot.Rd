\name{cover.plot}
\alias{cover.plot}
\alias{coverage}
\title{Unit-level coverage plots for beta parameters from MD EI model}

\description{Generates a plot of central credible intervals for the
  unit-level beta parameters from the Multinomial-Dirichlet ecological inference model
  (see \code{\link{ei.MD.bayes}}).} 

\usage{
cover.plot(object, row, column, x = NULL, CI = 0.95,
          medians = TRUE, col = NULL, ylim = c(0,1), 
          ylab, lty = par("lty"), lwd = par("lwd"), ...) 
}

\arguments{
  \item{object}{output from \code{\link{ei.MD.bayes}}}
  \item{row}{a character string specifying the row marginal of interest}
  \item{column}{a character string specifying the column marginal of
    interest}
  \item{x}{an optional covariate to index the units along the x-axis}
  \item{CI}{a fraction between 0 and 1 (defaults to 0.95), specifying
    the coverage of the central credible interval to be plotted for each unit}
  \item{medians}{a logical value specifying whether to plot the median
    (defaults to \code{TRUE}).  If \code{medians = FALSE}, the medians
    are not plotted.}
  \item{col}{an optional vector of colors to be passed to
    \code{plot} and \code{segments}.  If \code{col} is of
    length two, then the first color is used for \code{plot} and the
    second for \code{segments}.}
  \item{ylim}{an optional range for the y-axis (defaults to \code{c(0,1)}).}
  \item{ylab}{an optional label for the y-axis (defaults to
    \code{Proportion of row in column}).} 
  \item{lty}{an optional line type passed to \code{segments}.}
  \item{lwd}{an optional line width argument passed to
    \code{segments}.}
  \item{...}{additional arguments passed to \code{plot}.}
}
\value{
  A plot with vertical intervals indicating the central credible
  intervals for each ecological unit.  
}
\seealso{\code{plot}, \code{segments}, \code{par}}

\author{
  Olivia Lau <\email{olivia.lau@post.harvard.edu}>
}

\keyword{hplot}
