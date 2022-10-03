#' Tune in the peak information: peak position and peak scale
#'
#' Based on the identified peak position, more precise estimation of the peak
#' information, i.e., peak position and peak scale, can be got by this
#' function. The basic idea is to cut the segment of spectrum near the
#' identified peaks, and then do similar procedures as
#' [peakDetectionCWT()], but with more detailed scales around the
#' estimated peak scale.
#'
#' The majorPeakInfo or peakIndex and peakScale must be provided.
#'
#' @param ms the mass spectrometry spectrum
#' @param majorPeakInfo return of [identifyMajorPeaks()]
#' @param peakIndex the m/z index of the identified peaks
#' @param peakScale the scales of the identified peaks
#' @param maxScale the maximum scale allowed for the peak
#' @param \dots other parameters of used by [getLocalMaximumCWT()],
#' [getRidge()], [identifyMajorPeaks()]
#' @return \item{peakCenterIndex}{the updated peak center m/z index}
#' \item{peakScale}{the updated peak scale} \item{peakValue}{the corresponding
#' peak value}
#' @author Pan Du
#' @seealso [peakDetectionCWT()]
#' @references Du, P., Kibbe, W.A. and Lin, S.M. (2006) Improved peak detection
#' in mass spectrum by incorporating continuous wavelet transform-based pattern
#' matching, Bioinformatics, 22, 2059-2065.
#' @keywords methods
#' @export
#' @examples
#'
#'
#' data(exampleMS)
#' SNR.Th <- 3
#' peakInfo <- peakDetectionCWT(exampleMS, SNR.Th = SNR.Th)
#' majorPeakInfo <- peakInfo$majorPeakInfo
#' betterPeakInfo <- tuneInPeakInfo(exampleMS, majorPeakInfo)
#' plot(500:length(exampleMS), exampleMS[500:length(exampleMS)], type = "l", log = "x")
#' abline(v = betterPeakInfo$peakCenterIndex, col = "red")
#'
tuneInPeakInfo <- function(ms, majorPeakInfo = NULL, peakIndex = NULL, peakScale = NULL, maxScale = 128, ...) {
    if (!is.null(majorPeakInfo)) {
        if (!all(c("peakIndex", "peakCenterIndex", "peakScale") %in% names(majorPeakInfo))) stop("Format of majorPeakInfo is incorret!")
        peakIndex <- majorPeakInfo$peakIndex
        peakScale <- majorPeakInfo$peakScale[names(peakIndex)]
        peakCenterIndex <- majorPeakInfo$peakCenterIndex[names(peakIndex)]
        peakValue <- majorPeakInfo$peakValue[names(peakIndex)]
        peakSNR <- majorPeakInfo$peakSNR[names(peakIndex)]
    } else {
        if (is.null(peakIndex) | is.null(peakScale)) {
            stop("majorPeakInfo or peakIndex and peakScale should be provided!")
        }
        peakCenterIndex <- peakIndex
        peakSNR <- NULL
        peakValue <- NULL
    }

    peakName <- names(peakIndex)
    if (is.null(peakIndex)) peakName <- as.character(peakIndex)
    peakCenterIndex.new <- NULL
    peakScale.new <- NULL
    peakValue.new <- NULL
    unProcessedInd <- NULL
    for     (i in 1:length(peakIndex)) {
        peak.i <- peakIndex[i]
        # peak.i <- peakCenterIndex[i]
        peakScale.i <- peakScale[i]
        if (peakScale.i + 4 < maxScale) {
            scales.i <- seq(peakScale.i - 4, peakScale.i + 4, 0.5)
        } else {
            scales.i <- seq(maxScale - 10, maxScale, 0.5)
        }
        winSize.i <- 16 * max(scales.i)

        ## Select the start of the spectrum
        if (peak.i - winSize.i / 2 < 1) {
            winSize.i <- (peak.i - 1) * 2
            scales.i <- scales.i[scales.i < peak.i / 16]
            start.i <- 1
        } else {
            start.i <- peak.i - winSize.i / 2
        }

        ## Select the end of the spectrum
        if (peak.i + winSize.i / 2 > length(ms)) {
            winSize.i <- (length(ms) - peak.i) * 2
            scales.i <- scales.i[scales.i < (length(ms) - peak.i) / 16]
            end.i <- length(ms)
        } else {
            end.i <- peak.i + winSize.i / 2
        }
        if (length(scales.i) <= 1) {
            peakScale.new <- c(peakScale.new, peakScale[i])
            peakValue.new <- c(peakValue.new, peakValue[i])
            peakCenterIndex.new <- c(peakCenterIndex.new, peakCenterIndex[i])
            unProcessedInd <- c(unProcessedInd, i)
            next
        }

        ms.i <- ms[start.i:end.i]
        ## Perform Continuous Wavelet Transform
        wCoefs.i <- cwt(ms.i, scales = scales.i, wavelet = "mexh")

        ## -----------------------------------------
        ## Identify the local maximum by using a slide window
        localMax.i <- getLocalMaximumCWT(wCoefs.i, ...)
        colnames(localMax.i) <- colnames(wCoefs.i)

        ## -----------------------------------------
        ## Indentify the ridges from coarse level to more detailed levels
        ridgeList.i <- getRidge(localMax.i, gapTh = 3, skip = NULL, ...)

        ridgeName.i <- names(ridgeList.i)
        ridgeInfo.i <- matrix(as.numeric(unlist(strsplit(ridgeName.i, "_"))), nrow = 2)
        ridgeLevel.i <- ridgeInfo.i[1, ]
        newPeak.i <- ridgeInfo.i[2, ]

        # newPeak.i <- as.numeric(names(ridgeList.i))
        selInd.i <- which.min(abs(newPeak.i - winSize.i / 2))
        newRidgeLine.i <- ridgeList.i[[selInd.i]]
        newRidgeValue.i <- wCoefs.i[cbind(newRidgeLine.i, (1:length(newRidgeLine.i)) + ridgeLevel.i[selInd.i] - 1)]

        ## Find the maximum points of each route
        peakScaleInd.i <- which.max(newRidgeValue.i)
        newPeakValue.i <- max(newRidgeValue.i)
        peakScale.new <- c(peakScale.new, scales.i[peakScaleInd.i])
        peakValue.new <- c(peakValue.new, newPeakValue.i)
        peakCenterIndex.new <- c(peakCenterIndex.new, newRidgeLine.i[peakScaleInd.i] + start.i - 1)
    }
    if (!is.null(peakSNR)) {
        peakSNR.new <- peakSNR * peakValue.new / peakValue
    } else {
        peakSNR.new <- NULL
    }
    names(peakScale.new) <- names(peakValue.new) <- names(peakCenterIndex.new) <- peakName
    unProcessedPeak <- peakName[unProcessedInd]

    return(list(peakIndex = peakIndex, peakValue = peakValue.new, peakCenterIndex = peakCenterIndex.new, peakSNR = peakSNR.new, peakScale = peakScale.new, unProcessedPeak = unProcessedPeak))
}
