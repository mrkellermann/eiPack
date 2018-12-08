\name{read.betas}
\alias{read.betas}
\title{Function to read in eiMD parameter chains saved to disk}
\description{
  In \code{\link{ei.MD.bayes}}, users have the option to save parameter
  chains for the unit-level betas to disk rather than returning them to
  the workspace.  This function reconstructs the parameter chains by
  reading them back into R and producing either an array or an
  \code{mcmc} object.
}
\usage{
read.betas(rows, columns, units, dir = NULL, ret.mcmc = TRUE) 
}
\arguments{
  \item{rows}{a character vector of the row marginals to be read back
    in}
  \item{columns}{a character vector of the column marginals to be read
    back in}
  \item{units}{a character of numeric vector with the units to be read
    back in}
  \item{dir}{an optional character string identifying the directory in
    which parameter chains are stored (defaults to \code{getwd})}
  \item{ret.mcmc}{a logical value specifying whether to return the
    parameters as an \code{mcmc} object (defaults to \code{TRUE})}
}
\value{
  If \code{ret.mcmc = TRUE}, an \code{mcmc} object with row names
  corresponding to the parameter chains.  If \code{ret.mcmc = FALSE}, an
  array with dimensions named according to the selected \code{rows},
  \code{columns}, and \code{units}. 
}
\author{Olivia Lau <olivia.lau@post.harvard.edu>}
\seealso{\code{\link{ei.MD.bayes}},\code{mcmc}}
\keyword{IO}
\keyword{utilities}
