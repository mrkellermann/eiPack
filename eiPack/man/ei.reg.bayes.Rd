\name{ei.reg.bayes}
\alias{ei.reg.bayes}
\title{Ecological regression using Bayesian Normal regression}

\description{ Estimate an ecological regression using Bayesian normal 
regression. }

\usage{
ei.reg.bayes(formula, data, sample = 1000, weights = NULL, truncate=FALSE)
}

\arguments{
\item{formula}{An R formula object of the form \code{cbind(c1, c2, ...)
    ~ cbind(r1, r2, ...)}}
\item{data}{data frame containing the variables specified in formula}
\item{sample}{number of draws from the posterior}
\item{weights}{a vector of weights}
\item{truncate}{if TRUE, imposes a proper uniform prior on the unit
  hypercube for the coefficients; if FALSE, an improper uniform prior is assumed}
}

\value{ A list containing 
\item{call}{the call to \code{ei.reg.bayes}}
\item{draws}{A, \eqn{R \times C \times}{R x C x} sample array containing
  posterior draws for each population cell fraction}
}

\details{
  For \eqn{i \in 1,\ldots,C}{i in 1,...,C}, \eqn{C}{C} Bayesian regressions
  of the form \code{c_i ~ cbind(r1, r2, ...)} are
  performed.  See the documentation for \code{ei.reg} for the accounting
  identities and constancy assumption underlying this Bayesian linear
  model.

  The sampling density is given by \deqn{y|\beta, \sigma^2, X \sim
    N(X\beta, \sigma^2 I)}{y| beta, sigma^2, X ~ N(X beta, sigma^2*I)}

  The improper prior is \eqn{p(\beta,\sigma^2|X)\propto
    \sigma^{-2}}{p(beta,sigma^2|X) proportional to 1/sigma^2}.

  The proper prior is \eqn{p(\beta, \sigma^2|x) \propto I(\beta \in
    [0,1])\times \sigma^{-2}}{p(beta, sigma^2|x) proportional to I(beta
    in [0,1])* 1/sigma^2}.

}

\references{
  Leo Goodman. 1953. ``Ecological Regressions and the Behavior of Individuals.''
  \emph{American Sociological Review} 18:663--664.
}



\author{
  Olivia Lau <\email{olivia.lau@post.harvard.edu}> and Ryan T. Moore <\email{rtm@american.edu}>
}

\keyword{models}





