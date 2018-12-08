\name{mergeMD}
\alias{mergeMD}
\title{Combine output from multiple eiMD objects}

\description{Allows users to combine output from several chains
  output by \code{\link{ei.MD.bayes}}
}

\usage{
mergeMD(list, discard = 0) 
}

\arguments{
\item{list}{A list containing the names of multiple eiMD objects
  generated from the same model.}
\item{discard}{The number of draws to discard from the beginning of each chain.  
Default is to retain all draws.}
}
\value{
Returns an \code{eiMD} object of the same format as the input.
}


\references{
      Martyn Plummer, Nicky Best, Kate Cowles, and Karen Vines. 2002.
     \emph{Output Analysis and Diagnostics for MCMC (CODA)}. \url{
      https://CRAN.R-project.org/package=coda}.

	 Ori Rosen, Wenxin Jiang, Gary King, and Martin A. Tanner.
	 2001.  ``Bayesian and Frequentist Inference for Ecological
	 Inference: The \eqn{R \times C}{R x C} Case.'' \emph{Statistica
	 Neerlandica} 55: 134-156.  }



\author{
  Michael Kellermann <\email{mrkellermann@gmail.com}>
}

\seealso{\code{\link[eiPack]{ei.MD.bayes}}}

\keyword{utilities}





