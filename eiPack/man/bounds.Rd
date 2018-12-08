\name{bounds}
\alias{bounds}
\title{Deterministic bounds for units satisfying row thresholds}

\description{Calculates the deterministic bounds on the proportion of
  row members within a specified column.}

\usage{
bounds(formula, data, rows, column, excluded = NULL, 
    threshold = 0.9, total = NULL)
}

\arguments{
  \item{formula}{a formula of the form \code{cbind(col1, col2, ...) ~
      cbind(row1, row2, ...)}.  Column and row marginals must have the
    same total for each ecological unit.}
  \item{data}{a data frame containing the variables specified in
    \code{formula} and (optionally) \code{total}}
  \item{rows}{a character vector specifying the rows of interest}
  \item{column}{a character string specifying the column marginal of
    interest}
  \item{excluded}{an optional character string (or vector of character
    strings) specifying the columns to be excluded from the bounds
    calculation.  For example, if the quantity of interest is Democratic
    share of the two-party vote, non-voters would be excluded.}
  \item{threshold}{the minimum proportion of the unit that row members must
    comprise for the bounds to be calculated for the unit.  If
    \code{threshold = 0}, bounds will be calculated for all units.}
  \item{total}{if row and/or column marginals are given as proportions,
    \code{total} identifies the name of the variable in \code{data} containing the
    total number of individuals in each unit}
}
\value{ A list with elements 
  \item{bounds}{a list of deterministic bounds for all units in which
    row proportions meet the threshold}
  \item{intersection}{if the intersection of the deterministic bounding
    intervals is non-empty, the intersection is returned.  Otherwise,
    \code{NA} is returned.}
}
\seealso{\code{plot.bounds}}

\references{
  Otis Dudley Duncan and Beverley Davis.  1953.  ``An Alternative to Ecological
     Correlation.'' \emph{American Sociological Review} 18: 665-666.
}

\author{
  Ryan T. Moore <\email{rtm@american.edu}>
}

\keyword{models}
