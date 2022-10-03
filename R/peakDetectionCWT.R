#' The main function of peak detection by CWT based pattern matching
#'
#' This function is a wrapper of [cwt()],
#' [getLocalMaximumCWT()], [getRidge()],
#' [identifyMajorPeaks()]
#'
#'
#' @param ms the mass spectrometry spectrum
#' @param scales Scales of CWT. See [cwt()] for details.
#' @param SNR.Th SNR (Signal to Noise Ratio) threshold
#' @param nearbyPeak Determine whether to include the nearby small peaks of
#' major peaks. `TRUE` by default
#' @param peakScaleRange the scale range of the peak. larger than 5 by default.
#' @param amp.Th the minimum required relative amplitude of the peak (ratio to
#' the maximum of CWT coefficients)
#' @param minNoiseLevel the minimum noise level used in computing the SNR
#' @param ridgeLength the minimum highest scale of the peak in 2-D CWT
#' coefficient matrix
#' @param peakThr Minimal absolute intensity (above the baseline) of peaks to
#' be picked. If this value is provided, then the smoothing function
#' [signal::sgolayfilt()] will be called to estimate the local
#' intensity.(added based on the suggestion and code of Steffen Neumann)
#' @param tuneIn determine whether to tune in the parameter estimation of the
#' detected peaks
#' @param \dots other parameters used by [identifyMajorPeaks()].
#' Additionally, `fl` (filter length, with a default value of 1001) and
#' `forder` (filter order, with a default value of 2) are set and passed
#' to [signal::sgolayfilt()] when `peakThr` is given.
#' @return \item{majorPeakInfo}{ return of [identifyMajorPeaks()]}
#' \item{ridgeList}{return of [getRidge()]} \item{localMax}{ return
#' of [getLocalMaximumCWT()] } \item{wCoefs}{ 2-D CWT coefficient
#' matrix, see [cwt()] for details.}
#' @author Pan Du, Simon Lin
#' @seealso [cwt()], [getLocalMaximumCWT()],
#' [getRidge()], [identifyMajorPeaks()]
#' @references Du, P., Kibbe, W.A. and Lin, S.M. (2006) Improved peak detection
#' in mass spectrum by incorporating continuous wavelet transform-based pattern
#' matching, Bioinformatics, 22, 2059-2065.
#' @keywords methods
#' @export
#' @examples
#'
#' data(exampleMS)
#' SNR.Th <- 3
#' peakInfo <- peakDetectionCWT(exampleMS, SNR.Th = SNR.Th)
#' majorPeakInfo <- peakInfo$majorPeakInfo
#' peakIndex <- majorPeakInfo$peakIndex
#' plotPeak(exampleMS, peakIndex, main = paste("Identified peaks with SNR >", SNR.Th))
#'
#' ## In some cases, users may want to add peak filtering based on the absolute peak amplitude
#' peakInfo <- peakDetectionCWT(exampleMS, SNR.Th = SNR.Th, peakThr = 500)
#' majorPeakInfo <- peakInfo$majorPeakInfo
#' peakIndex <- majorPeakInfo$peakIndex
#' plotPeak(exampleMS, peakIndex, main = paste("Identified peaks with SNR >", SNR.Th))
#'
peakDetectionCWT <- function(ms, scales = c(1, seq(2, 30, 2), seq(32, 64, 4)), SNR.Th = 3, nearbyPeak = TRUE, peakScaleRange = 5, amp.Th = 0.01, minNoiseLevel = amp.Th / SNR.Th, ridgeLength = 24, peakThr = NULL, tuneIn = FALSE, ...) {
    otherPar <- list(...)
    if (minNoiseLevel > 1) names(minNoiseLevel) <- "fixed"
    ## Perform Continuous Wavelet Transform
    wCoefs <- cwt(ms, scales = scales, wavelet = "mexh")

    ## Attach the raw data as the zero level of decomposition
    wCoefs <- cbind(as.vector(ms), wCoefs)
    colnames(wCoefs) <- c(0, scales)

    ## -----------------------------------------
    ## Identify the local maximum by using a slide window
    ## The size of slide window changes over different levels, with the coarse level have bigger window size
    if (is.null(amp.Th)) {
        amp.Th <- 0
    } else {
        if (is.null(names(amp.Th))) {
            amp.Th <- max(wCoefs) * amp.Th
        } else if (names(amp.Th) != "fixed") {
            amp.Th <- max(wCoefs) * amp.Th
        }
    }
    localMax <- getLocalMaximumCWT(wCoefs, amp.Th = amp.Th)
    colnames(localMax) <- colnames(wCoefs)

    ## In order to fastern the calculation, we can filter some local maxima with small amplitude
    ## In this case a baseline estimation was performed.
    if (!is.null(peakThr)) {
        if (!requireNamespace("signal", quietly = TRUE)) {
            stop('The use of peakThr= argument in MassSpecWavelet::peakDetectionCWT() requires to install the "signal" package. Please use install.packages("signal")')
        }
        if ("fl" %in% names(otherPar)) {
            filterLength <- otherPar$fl
        } else {
            filterLength <- 1001
        }
        if (filterLength %% 2 == 0) {
            warning("filter length in peakDetectionCWT(fl=) needs to be odd (increasing it by 1)")
            filterLength <- filterLength + 1
        }
        if ("forder" %in% names(otherPar)) {
            fOrder <- otherPar$forder
        } else {
            fOrder <- 2
        }
        ## Baseline estimation using Savitzky Golay Filter
        ## this part was added by Steffen Neumann
        sg <- signal::sgolayfilt(ms, p = fOrder, n = filterLength)
        localMax[(ms - sg) < peakThr, ] <- 0
    }
    ## remove the parameters in otherPar that were passed to the Savitzky-Golay filter
    otherPar <- otherPar[!(names(otherPar) %in% c("fl", "forder", "dorder"))]

    ## -----------------------------------------
    ## Indentify the ridges from coarse level to more detailed levels
    ridgeList <- getRidge(localMax, gapTh = 3, skip = 2)

    ## -----------------------------------------
    ## Indentify the major peaks and their nearby peaks
    majorPeakInfo <- do.call(identifyMajorPeaks, c(list(
        ms = ms, ridgeList = ridgeList, wCoefs = wCoefs, SNR.Th = SNR.Th, peakScaleRange = peakScaleRange,
        nearbyPeak = nearbyPeak, minNoiseLevel = minNoiseLevel, ridgeLength = ridgeLength
    ), otherPar))

    if (tuneIn) {
        refinedPeakInfo <- tuneInPeakInfo(ms, majorPeakInfo)
        return(list(majorPeakInfo = refinedPeakInfo, ridgeList = ridgeList, localMax = localMax, wCoefs = wCoefs[, -1], oldPeakInfo = majorPeakInfo))
    } else {
        return(list(majorPeakInfo = majorPeakInfo, ridgeList = ridgeList, localMax = localMax, wCoefs = wCoefs[, -1]))
    }
}
