\name{peakDetectionCWT}
\alias{peakDetectionCWT}
\title{ The main function of peak detection by CWT based pattern matching }
\description{
  This function is a wrapper of \code{\link{cwt}}, \code{\link{getLocalMaximumCWT}}, \code{\link{getRidge}}, \code{\link{identifyMajorPeaks}}
}
\usage{
peakDetectionCWT(
 ms,
 scales = c(1, seq(2, 30, 2), seq(32, 64, 4)),
 SNR.Th = 3,
 nearbyPeak = TRUE,
 peakScaleRange = 5,
 amp.Th = 0.01,
 minNoiseLevel = amp.Th/SNR.Th,
 ridgeLength = 24,
 peakThr=NULL,
 tuneIn = FALSE,
 ...
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ms}{ the mass spectrometry spectrum }
  \item{scales}{ Scales of CWT. See \code{\link{cwt}} for details. }
  \item{SNR.Th}{ SNR (Signal to Noise Ratio) threshold}
  \item{nearbyPeak}{ Determine whether to include the nearby small peaks of major peaks. \code{TRUE} by default}
  \item{peakScaleRange}{ the scale range of the peak. larger than 5 by default. }
  \item{amp.Th}{ the minimum required relative amplitude of the peak (ratio to the maximum of CWT coefficients) }
  \item{minNoiseLevel}{ the minimum noise level used in computing the SNR }
  \item{ridgeLength}{ the minimum highest scale of the peak in 2-D CWT coefficient matrix }
  \item{peakThr}{ Minimal absolute intensity (above the baseline) of peaks to be picked. If this value is provided, then the smoothing function \code{\link[signal]{sgolayfilt}} will be called to estimate the local intensity.(added based on the suggestion and code of Steffen Neumann)}
  \item{tuneIn}{ determine whether to tune in the parameter estimation of the detected peaks }
  \item{\dots}{other parameters used by \code{\link{identifyMajorPeaks}}. Additionally, \code{fl} (filter length, with a default value of 1001) and \code{forder} (filter order, with a default value of 2) are set and passed to \code{\link[signal]{sgolayfilt}} when \code{peakThr} is given.}
}

\value{
  \item{majorPeakInfo}{ return of \code{\link{identifyMajorPeaks}}}
  \item{ridgeList}{return of \code{\link{getRidge}}}
  \item{localMax}{ return of \code{\link{getLocalMaximumCWT}} }
  \item{wCoefs}{ 2-D CWT coefficient matrix, see \code{\link{cwt}} for details.}
}

\references{ 
Du, P., Kibbe, W.A. and Lin, S.M. (2006) Improved peak detection in mass spectrum by incorporating continuous wavelet transform-based pattern matching, Bioinformatics, 22, 2059-2065.
}
\author{Pan Du, Simon Lin}

\seealso{ \code{\link{cwt}}, \code{\link{getLocalMaximumCWT}}, \code{\link{getRidge}}, \code{\link{identifyMajorPeaks}} }
\examples{
	data(exampleMS)
	SNR.Th <- 3
	peakInfo <- peakDetectionCWT(exampleMS, SNR.Th=SNR.Th)
	majorPeakInfo = peakInfo$majorPeakInfo
	peakIndex <- majorPeakInfo$peakIndex
	plotPeak(exampleMS, peakIndex, main=paste('Identified peaks with SNR >', SNR.Th)) 
	
	## In some cases, users may want to add peak filtering based on the absolute peak amplitude
	peakInfo <- peakDetectionCWT(exampleMS, SNR.Th=SNR.Th, peakThr=500)
	majorPeakInfo = peakInfo$majorPeakInfo
	peakIndex <- majorPeakInfo$peakIndex
	plotPeak(exampleMS, peakIndex, main=paste('Identified peaks with SNR >', SNR.Th)) 	
}
\keyword{methods}

