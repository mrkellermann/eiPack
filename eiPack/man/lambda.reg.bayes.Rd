\name{lambda.reg.bayes}
\alias{lambda.reg.bayes}
\title{Calculate shares using data from Bayesian regression model}

\description{Calculates the population share of row members in selected columns}

\usage{
lambda.reg.bayes(object, columns, ret.mcmc = TRUE) 
}

\arguments{
\item{object}{An R object of class \code{eiRegBayes}, the output from \code{\link{ei.reg.bayes}}}
\item{columns}{a character vector indicating which column marginals to be included in calculating the shares}
\item{ret.mcmc}{If TRUE, posterior shares are returned as an \code{mcmc}
  object.}
}
\value{
  If \code{ret.mcmc = TRUE}, draws are returned as an \code{mcmc} object
  with dimensions sample \eqn{\times C}{x C}.  If \code{ret.mcmc =
    FALSE}, draws are returned as an array with dimensions \eqn{R
    \times C \times }{R x C x} samples array.  
}

\seealso{\code{\link{ei.reg.bayes}}}

\author{
  Ryan T. Moore <\email{rtm@american.edu}>
}

\keyword{models}
