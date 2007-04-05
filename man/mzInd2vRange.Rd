\name{mzInd2vRange}
\alias{mzInd2vRange}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Match m/z index to m/z value with a certain error range }
\description{
Match m/z index to m/z value with a certain error range 
}

\usage{
mzInd2vRange(mzInd, error = 0.003)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mzInd}{ a vector of m/z index }
  \item{error}{ error range }
}

\value{
 return a vector of sorted m/z values
}

\author{Pan Du}

\seealso{ \code{\link{mzV2indRange}}}

\keyword{methods}
