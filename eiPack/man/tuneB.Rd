\name{tuneB}
\alias{tuneB}

\title{Tuning parameters for the precinct level parameters in the RxC EI 
model}

\description{ A vector containing tuning parameters for the precinct level 
parameters in the RxC EI model.  
}

\usage{data(tuneB)}

\format{ A vector of length 3 x 2 x 150 containing the precinct level tuning 
parameters for the redistricting sample data.  
}

\examples{
data(tuneB)
tuneB <- array(tuneB[[1]], dim = c(3, 2, 150))
}
\keyword{datasets}
