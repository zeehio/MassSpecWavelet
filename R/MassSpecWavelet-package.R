#' Peak detection of mass spectrum by Wavelet transform based methods
#'
#' Process Mass Spectrum (MS) by Wavelet Transforms-based algorithms
#'
#' \tabular{ll}{ Package: \tab MassSpecWavelet\cr Type: \tab Package\cr
#' Version: \tab 1.0.4\cr Date: \tab 2007-04-05\cr License: \tab GPL 2 or
#' newer\cr } MassSpecWavelet R package is aimed to process Mass Spectrometry
#' (MS) data mainly based on Wavelet Transforms. The current version only
#' supports the peak detection based on Continuous Wavelet Transform (CWT).
#' Future versions will include more functions covering entire MS data
#' processes.
#'
#' @author Pan Du, Simon Lin
#'
#' Maintainer: Pan Du <dupan@@northwestern.edu>
#' @references Du, P., Kibbe, W.A. and Lin, S.M. (2006) Improved peak detection
#' in mass spectrum by incorporating continuous wavelet transform-based pattern
#' matching, Bioinformatics, 22, 2059-2065.
#' @keywords package
#' @examples
#'
#' data(exampleMS)
#' SNR.Th <- 3
#' peakInfo <- peakDetectionCWT(exampleMS, SNR.Th = SNR.Th)
#' majorPeakInfo <- peakInfo$majorPeakInfo
#' peakIndex <- majorPeakInfo$peakIndex
#' plotPeak(exampleMS, peakIndex, main = paste("Identified peaks with SNR >", SNR.Th))
#'
"_PACKAGE"


#' An example mass spectrum
#'
#' An example mass spectrum from CAMDA 2006. All-in-1 Protein Standard II
#' (Ciphergen Cat. \# C100-0007) were measured on Ciphergen NP20 chips. There
#' are 7 polypeptides in the sample with m/z values of 7034, 12230, 16951,
#' 29023, 46671, 66433, 147300.
#'
#'
#' @name exampleMS
#' @docType data
#' @format A numeric vector represents the mass spectrum with equal sample
#' intervals.
#' @source CAMDA, CAMDA 2006 Competition Data Set. 2006, http://camda.duke.edu.
#' @keywords datasets
NULL
