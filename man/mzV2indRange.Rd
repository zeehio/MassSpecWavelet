\name{mzV2indRange}
\alias{mzV2indRange}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Match m/z value to m/z index with a certain error range}
\description{
Match m/z value to m/z index with a certain error range
}
\usage{
mzV2indRange(mzV, error = 0.003)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mzV}{ a vector of m/z value }
  \item{error}{ error range }
}
\value{
return a vector of sorted m/z indexes
}
\author{Pan Du}

\seealso{ \code{\link{mzInd2vRange}}}

\keyword{methods}
