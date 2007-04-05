\name{MassSpecWavelet-package}
\alias{MassSpecWavelet-package}
\alias{MassSpecWavelet}
\docType{package}
\title{Peak detection of mass spectrum by Wavelet transform based methods}
\description{
Process Mass Spectrum (MS) by Wavelet Transforms-based algorithms
}
\details{
\tabular{ll}{
Package: \tab MassSpecWavelet\cr
Type: \tab Package\cr
Version: \tab 1.0.2\cr
Date: \tab 2007-03-02\cr
License: \tab GPL 2 or newer\cr
}
MassSpecWavelet R package is aimed to process Mass Spectrometry (MS) data 
mainly based on Wavelet Transforms. The current version only supports the 
peak detection based on Continuous Wavelet Transform (CWT). Future versions 
will include more functions covering entire MS data processes.

}
\author{
Pan Du, Simon Lin

Maintainer: Pan Du <dupan@northwestern.edu>
}
\references{
Pan Du, Warren A. Kibbe, Simon M. Lin, Improved peak detection in mass spectrum by incorporating continuous wavelet transform-based pattern matching, Bioinformatics Advance Access published online on July 4, 2006
(http://bioinformatics.oxfordjournals.org/cgi/content/abstract/btl355?ijkey=HfCb7kHvQBY6HE6&keytype=ref)
}
\keyword{ package }

\examples{
data(exampleMS)
SNR.Th <- 3
peakInfo <- peakDetectionCWT(exampleMS, SNR.Th=SNR.Th)
majorPeakInfo = peakInfo$majorPeakInfo
peakIndex <- majorPeakInfo$peakIndex
plotPeak(exampleMS, peakIndex, main=paste('Identified peaks with SNR >', SNR.Th)) 
}
