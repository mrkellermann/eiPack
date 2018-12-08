\name{tuneMD}
\alias{tuneMD}
\title{Generate tuning parameters for MD model}

\description{An adaptive algorithm to generate tuning parameters for the MCMC 
algorithm implemented in \code{\link{ei.MD.bayes}}.  Since we are 
drawing each parameter one at a time, target acceptance rates are between 0.4 to 0.6.}

\usage{
tuneMD(formula, covariate = NULL, data, ntunes = 10, 
    totaldraws = 10000, ...)
}

\arguments{
\item{formula}{A formula of the form \code{cbind(col1, col2, ...) ~ cbind(row1, row2, ...)} with rows as the predictor and columns 
as the response}
\item{covariate}{An R formula for the optional covariate in the form
  \code{~ x}}
\item{data}{data frame containing the variables specified in \code{formula} and
\code{covariate}}
\item{ntunes}{number of times to iterate the tuning algorithm}
\item{totaldraws}{number of iterations for each tuning run}
\item{...}{additional arguments passed to \code{\link{ei.MD.bayes}}}
}

\value{
A list containing matrices of tuning parameters.
}

\seealso{\code{\link{ei.MD.bayes}}}

\author{
  Olivia Lau <\email{olivia.lau@post.harvard.edu}>
}

\keyword{iteration}
\keyword{utilities}





