\name{plotLocalMax}
\alias{plotLocalMax}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot the local maximum matrix}
\description{
Plot the local maximum matrix of 2-D CWT coefficients returned by \code{\link{getLocalMaximumCWT}}
}
\usage{
plotLocalMax(
 localMax,
 wCoefs = NULL,
 range = c(1, nrow(localMax)),
 colorMap = "RYB",
 main = NULL, cex = 3, pch = ".", ...
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{localMax}{local maximum matrix of 2-D CWT coefficients returned by \code{\link{getLocalMaximumCWT}}}
  \item{wCoefs}{ 2-D CWT coefficients }
  \item{range}{ plot range of m/z index }
  \item{colorMap}{ the colormap used in plotting the points }
  \item{main}{ parameter of \code{\link{plot}} }
  \item{cex}{ parameter of \code{\link{plot}} }
  \item{pch}{ parameter of \code{\link{plot}} }
  \item{\dots}{ other parameters of \code{\link{points}} }
}

\author{Pan Du}

\seealso{ \code{\link{getLocalMaximumCWT}} }
\examples{
	data(exampleMS)
	scales <- seq(1, 64, 3)
	wCoefs <- cwt(exampleMS[5000:11000], scales=scales, wavelet='mexh')
	
	localMax <- getLocalMaximumCWT(wCoefs)
	plotLocalMax(localMax)
}
\keyword{hplot}
