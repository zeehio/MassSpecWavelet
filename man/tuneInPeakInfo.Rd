\name{tuneInPeakInfo}
\alias{tuneInPeakInfo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Tune in the peak information: peak position and peak scale }
\description{
Based on the identified peak position, more precise estimation of the peak information, i.e., peak position and peak scale, can be got by this function. The basic idea is to cut the segment of spectrum near the identified peaks, and then do similar procedures as \code{\link{peakDetectionCWT}}, but with more detailed scales around the estimated peak scale.
}
\usage{
tuneInPeakInfo(ms, majorPeakInfo = NULL, peakIndex = NULL,
               peakScale = NULL, maxScale = 128, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ms}{ the mass spectrometry spectrum }
  \item{majorPeakInfo}{ return of \code{\link{identifyMajorPeaks}} }
  \item{peakIndex}{ the m/z index of the identified peaks}
  \item{peakScale}{ the scales of the identified peaks }
  \item{maxScale}{ the maximum scale allowed for the peak }
  \item{\dots}{ other parameters of used by \code{\link{getLocalMaximumCWT}}, \code{\link{getRidge}}, \code{\link{identifyMajorPeaks}} }
}
\details{
 The majorPeakInfo or peakIndex and peakScale must be provided.
}
\value{
  \item{peakCenterIndex}{the updated peak center m/z index}
  \item{peakScale}{the updated peak scale}
  \item{peakValue}{the corresponding peak value}
}
\references{ 
Du, P., Kibbe, W.A. and Lin, S.M. (2006) Improved peak detection in mass spectrum by incorporating continuous wavelet transform-based pattern matching, Bioinformatics, 22, 2059-2065.
}
\author{ Pan Du }

\seealso{ \code{\link{peakDetectionCWT}} }
\examples{

	data(exampleMS)
	SNR.Th <- 3
	peakInfo <- peakDetectionCWT(exampleMS, SNR.Th=SNR.Th)
	majorPeakInfo <- peakInfo$majorPeakInfo
	betterPeakInfo <- tuneInPeakInfo(exampleMS, majorPeakInfo)
	plot(500:length(exampleMS), exampleMS[500:length(exampleMS)], type='l', log='x')
	abline(v=betterPeakInfo$peakCenterIndex, col='red')

}
\keyword{methods}
