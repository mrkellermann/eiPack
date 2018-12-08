\name{ei.reg}
\alias{ei.reg}
\title{Ecological regression}

\description{ Estimate an ecological regression using least squares. }

\usage{
ei.reg(formula, data, ...)
}

\arguments{
\item{formula}{An R formula object of the form \code{cbind(c1, c2, ...)
    ~ cbind(r1, r2, ...)}} 
\item{data}{data frame containing the variables specified in \code{formula}}
\item{\dots}{Additional arguments passed to \code{\link[stats]{lm}}. }  
}


\value{A list containing 
    \item{call}{the call to \code{ei.reg}}
    \item{coefficients}{an \eqn{R \times C}{R x C} matrix of estimated
      population cell fractions}
    \item{se}{an \eqn{R \times C}{R x C} matrix of standard errors
    for \code{coefficients}.}
    \item{cov.matrices}{A list of the \eqn{C}{C} scaled variance-covariance
      matrices for each of the ecological regressions}
}

\details{
  For \eqn{i \in 1,\ldots,C}{i in 1,...,C}, C regressions of the form
  \code{c_i ~ cbind(r1, r2, ...)} are performed.

  These regressions make use of the accounting identities
  and the constancy assumption, that \eqn{\beta_{rci} =
    \beta_{rc}}{beta_rci = beta_rc} for all \eqn{i}{i}.

  The accounting identities include
  \itemize{
    \item{--}{defining the population cell fractions
    \eqn{\beta_{rc}}{beta_rc} such that \eqn{\sum_{c=1}^{C} \beta_{rc}
      =  1}{sum_{c=1}^{C} beta_rc  =  1} for every \eqn{r}{r}}
  \item{--}{\eqn{\sum_{c=1}^{C} \beta_{rci}  =  1}{sum_{c=1}^{C} beta_rci  =
      1} for \eqn{r = 1, \ldots, R}{r = 1,...,R} and \eqn{i = 1, \ldots,
      n}{i = 1,...,n}}
  \item{--}{\eqn{T_{ci} = \sum_{r=1}^R \beta_{rci}X_{ri}}{T_ci = sum_{r=1}^R
      beta_rci X_ri} for \eqn{c = 1,\ldots,C}{c = 1,...,C} and \eqn{i =
      1\ldots,n}{i = 1,...,n}}
  }

  Then regressing \deqn{T_{ci} = \beta_{rc} X_{ri} + \epsilon_{ci}}{T_ci
    = beta_rc X_ri + epsilon_ci} for \eqn{c = 1,\dots,C}{c = 1,...C}
  recovers the population parameters \eqn{\beta_{rc}}{beta_rc} when the
  standard linear regression assumptions apply, including
  \eqn{E[\epsilon_{ci}] = 0}{E[epsilon_ci] = 0} and
  \eqn{Var[\epsilon_{ci}] = \sigma_c^2}{Var[epsilon_ci] = sigma_c^2} for
  all \eqn{i}{i}.
}

\references{
   Leo Goodman. 1953. ``Ecological Regressions and the Behavior of Individuals.''
   \emph{American Sociological Review} 18:663--664.
 }

\author{
  Olivia Lau <\email{olivia.lau@post.harvard.edu}> and Ryan T. Moore <\email{rtm@american.edu}>
}

\keyword{models}





