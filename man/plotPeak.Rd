\name{plotPeak}
\alias{plotPeak}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot the identified peaks over the spectrum}
\description{
Plot the identified peaks over the spectrum. The identified peaks are returned by \code{\link{peakDetectionCWT}} or \code{\link{identifyMajorPeaks}}
}
\usage{
plotPeak(ms, peakIndex = NULL, mz = 1:length(ms),
         range = c(min(mz), max(mz)), method = c("p", "l"), main = NULL, log = "", ...
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ms}{ the MS spectrum }
  \item{peakIndex}{ m/z indexes of the identified peaks }
  \item{mz}{ m/z value correspond to m/z index }
  \item{range}{ the plot range of m/z value}
  \item{method}{ plot method of the identified peaks. method 'p' plot circles on the peaks; method 'l' add vertical lines over the peaks.  }
  \item{main}{ parameter of \code{\link{plot}}}
  \item{log}{ parameter of \code{\link{plot}} }
  \item{\dots}{ other parameters of \code{\link{points}}}
}

\author{Pan Du}
\value{
No value is returned; this function is called for its side effects (plot).
}

\seealso{\code{\link{peakDetectionCWT}}, \code{\link{identifyMajorPeaks}}}
\examples{
	data(exampleMS)
	SNR.Th <- 3
	peakInfo <- peakDetectionCWT(exampleMS, SNR.Th=SNR.Th)
	majorPeakInfo = peakInfo$majorPeakInfo
	peakIndex <- majorPeakInfo$peakIndex
	plotPeak(exampleMS, peakIndex, main=paste('Identified peaks with SNR >', SNR.Th)) 

}
\keyword{hplot}
