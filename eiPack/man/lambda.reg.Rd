\name{lambda.reg}
\alias{lambda.reg}
\title{Calculate shares using data from regression model}

\description{Calculates the population share of row members in a particular column}

\usage{
lambda.reg(object, columns) 
}

\arguments{
\item{object}{An R object of class \code{eiReg}, the output from \code{\link{ei.reg}}}
\item{columns}{a character vector of column names to be included in calculating
the shares}
}
\value{
Returns a list with the following elements 
\item{call}{the call to \code{lambda.reg}}
\item{lambda}{an \eqn{R \times k}{R x k} matrix 
where \eqn{k}{k} is the
number of columns included in the share calculation}
\item{se}{standard errors calculated using the delta method as implemented 
in the library \code{msm}}
}

\details{
Standard errors are calculated using the delta method as implemented in
the library \code{msm}.  The arguments passed to
\code{deltamethod} in \code{msm} include 
  \itemize{
    \item{\code{g}}{a list of transformations of the form \code{~ x1 / (x1 + x2 +
      + ... + xk)}, \code{~ x2 / (x1 + x2 + ... + xk)}, etc.}. Each
  \eqn{x_c}{xc} is the estimated proportion of all row members in column
  \eqn{c}{c}, \eqn{\hat{\beta}_{rc}}{beta_rc}
    \item{\code{mean}}{the estimated proportions of the row members in the
      specified columns, as a proportion of the total number of row
      members, \eqn{(\hat{\beta}_{r1}, \hat{\beta}_{r2}, ...,
	\hat{\beta}_{rk})}{(beta_r1, beta_r2, ..., beta_rk)}.} 
    \item{\code{cov}}{a diagonal matrix with the estimated variance of each
      \eqn{\hat{\beta}_{rc}}{beta_rc} on the diagonal.  Each column
      marginal is assumed to be independent, such that the off-diagonal
      elements of this matrix are zero.  Estimates come from
      \code{object$cov.matrices}, the estimated covariance matrix from
      the regression of the relevant column.  Thus,
    }
  }
      \tabular{cccccc}{
      cov \tab = \tab \eqn{Var(\hat{\beta}_{r1})}{Var(beta_r1)} \tab 0
      \tab 0 \tab \eqn{\ldots}{...} \cr 
\tab \tab 0 \tab \eqn{Var(\hat{\beta}_{r2})}{Var(beta_r2)} \tab 0 \tab \eqn{\ldots}{...} \cr
\tab \tab 0 \tab 0 \tab \eqn{Var(\hat{\beta}_{r3})}{Var(beta_{r3})} \tab \eqn{\ldots}{...} \cr
\tab \tab \eqn{\vdots}{...} \tab \eqn{\vdots}{...} \tab \eqn{\vdots}{...} \tab \eqn{\ddots}{...}\cr
  }
}

\seealso{\code{\link{ei.reg}}}

\author{
  Ryan T. Moore <\email{rtm@american.edu}>
}

\keyword{models}
