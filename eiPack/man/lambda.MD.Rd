\name{lambda.MD}
\alias{lambda.MD}
\title{Calculate shares using data from MD model}

\description{Calculates the population share of row members in a
  particular column as a proportion of the total number of row members in the selected subset of columns.}

\usage{
lambda.MD(object, columns, ret.mcmc = TRUE)
}

\arguments{
\item{object}{an R object of class \code{eiMD}, output from \code{\link{ei.MD.bayes}}}
\item{columns}{a character vector of column names to be included in
  calculating the shares} 
\item{ret.mcmc}{a logical value indicating how the samples from the posterior 
should be returned.  If \code{TRUE} (default), samples are returned as \code{mcmc} 
objects.  If \code{FALSE}, samples are returned as arrays. }
}
\value{
Returns either a ((\eqn{R}{R} * included columns) \eqn{\times}{x}
samples) matrix as an \code{mcmc} object or a (\eqn{R \times}{R x}
included columns \eqn{\times}{x} samples) array. 
}

\details{This function allows users to define subpopulations within the
  data and calculate the proportion of individuals within each of the
  columns that defines that subpopulation.  For example, if the model
  includes the groups Democrat, Republican, and Unaffiliated, the
  argument \code{columns = c(``Democrat", ``Republican")} will calculate
  the two-party shares of Democrats and Republicans for each row.
}

\seealso{\code{\link{ei.MD.bayes}}}

\author{
  Michael Kellermann <\email{mrkellermann@gmail.com}> and Olivia Lau 
<\email{olivia.lau@post.harvard.edu}>
}

\keyword{models}





