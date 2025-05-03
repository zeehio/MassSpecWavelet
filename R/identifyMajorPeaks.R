#' Identify peaks based on the ridges in 2-D CWT coefficient matrix
#'
#' Identify the peaks based on the ridge list (returned by
#' [getRidge()]) in 2-D CWT coefficient matrix and estimated Signal
#' to Noise Ratio (SNR). The criteria for peak identification is described in
#' the Details section.
#'
#' The ridge list may return peaks than have to be filtered out. This function
#' filters the peaks according to the following rules. All rules must pass for a
#' peak to be identified as such.
#' 
#' \itemize{
#'   \item{The maximum scale of the peak ridge should be larger than `ridgeLength`. 
#'   If `nearbyPeak=TRUE`, all peaks at less than `nearbyWinSize` points from a
#'   peak that fullfills this rule are also considered.}
#'   \item{The SNR of the peak must be larger than `SNR.Th`.}
#'   \item{The peak should not appear at the first or last `excludeBoundariesSize` points.}
#' }
#' 
#' To debug and diagnose why good peaks get filtered, you may want to set `ridgeLength=0`,
#' `SNR.Th=0` and/or `excludeBoundariesSize=0` to disable each of the filtering
#' criteria.
#' 
#' ## SNR estimation
#' 
#' The SNR is defined as \eqn{SNR = \frac{signal}{noise}}{SNR=signal/noise}. Both
#' signal and noise values need to be estimated for each peak.
#' 
#' The "signal" is estimated as the maximum wavelet coefficient obtained in the
#' corresponding peak ridge, considering all the scales within `peakScaleRange`.
#' 
#' The "noise" is estimated differently depending on the `SNR.method`. All methods
#' use a window of data points of size `2 * winSize.noise + 1` centered at the peak
#' to make the noise estimation. Here is how the noise is estimated depending on
#' the `SNR.method` value:
#' 
#' \itemize{
#'   \item{`"quantile"`: The "noise" is the 95% quantile of the absolute value of the wavelet 
#'   coefficients at scale 1 in the window.}
#'   \item{`"sd"`: The "noise" is the standard deviation of the absolute value of the wavelet
#'   coefficients at scale 1 in the window.}
#'   \item{`"mad"`: The "noise" is the [mad()] with `center=0` of the absolute value of the
#'   wavelet coefficients at scale 1 in the window.}
#'   \item{`"data.mean"`: The "noise" is the mean value of the ms spectrum in the window.}
#'   \item{`"data.mean.quant"`: The "noise" is the mean value of the ms spectrum in the window,
#'   but only considering values below the 95% quantile in the window.}
#' }
#' 
#' 
#' If the obtained noise estimation is below the minimum noise level, that minimum
#' is used as the noise estimation instead. Check `minNoiseLevel` for further details
#' on how the minimum noise level is defined.
#' 
#' Using the estimated "signal" and "noise", we compute the `peakSNR` value for each peak.
#' 
#'
#' @param ms the mass spectrometry spectrum
#' @param ridgeList returned by [getRidge()]
#' @param wCoefs 2-D CWT coefficients as obtained by [cwt()].
#' @param scales scales of CWT, by default it is the colnames of wCoefs
#' @param SNR.Th threshold of SNR
#' @param peakScaleRange the CWT scale range of the peak, used to estimate the
#' signal of the SNR. See Details. If a single value is given then
#' all scales larger than the value will be considered. If two values are given
#' only the scales between those values will be considered.
#' @param ridgeLength the maximum ridge scale of the major peaks.
#' @param nearbyPeak determine whether to include the small peaks close to
#' large major peaks. See Details.
#' @param nearbyWinSize the window size to determine the nearby peaks. Only
#' effective when `nearbyPeak=TRUE`.
#' @param winSize.noise the local window size to estimate the noise level.
#' @param SNR.method method to estimate the noise level. See Details.
#' @param minNoiseLevel the minimum noise level used in calculating SNR.
#' This value should be zero or positive. If the number is smaller than one, it
#' is assumed to be a fraction of the largest wavelet coefficient in the data.
#' Otherwise it is assumed to be the actual noise level. If you want to fix the
#' actual noise level to a value smaller than one, you should name the value as fixed
#' as in `minNoiseLevel = c("fixed"= 0.5)`. See details.
#' @param excludeBoundariesSize number of points at each boundary of the ms
#' signal that will be excluded in search for peaks to avoid boundary effects.
#' 
#' @return Return a list with following elements:
#' 
#' \describe{
#' \item{peakIndex}{the m/z indexes of the identified peaks}
#' \item{peakCenterIndex}{the m/z indexes of peak centers, which correspond to
#' the maximum on the ridge. `peakCenterIndex` includes all the peaks, not just
#' the identified major peaks.}
#' \item{peakValue}{the CWT coefficients (the maximum on the ridge)
#' corresponding to peakCenterIndex}
#' \item{peakSNR}{the SNR of the peak, which is the ratio of peakValue and noise
#' level}
#' \item{peakScale}{the estimated scale of the peak, which corresponds to the `peakCenterIndex`}
#' \item{potentialPeakIndex}{the m/z indexes of all potential peaks, which
#' satisfy all requirements of a peak without considering its SNR. Useful, if
#' you want to change to a lower SNR threshold later.}
#' \item{allPeakIndex}{the m/z indexes of all the peaks, whose order is the 
#' same as `peakCenterIndex`, `peakCenterValue`, `peakSNR` and `peakScale`.}
#' }
#' 
#' All of these return elements have peak names, which are the same as the
#' corresponding peak ridges. see [getRidge()] for details.
#' 
#' @author Pan Du, Simon Lin
#' @seealso [peakDetectionCWT()], [tuneInPeakInfo()]
#' @references Du, P., Kibbe, W.A. and Lin, S.M. (2006) Improved peak detection
#' in mass spectrum by incorporating continuous wavelet transform-based pattern
#' matching, Bioinformatics, 22, 2059-2065.
#' @keywords methods
#' @export
#' @examples
#'
#' data(exampleMS)
#' scales <- seq(1, 64, 3)
#' wCoefs <- cwt(exampleMS, scales = scales, wavelet = "mexh")
#'
#' localMax <- getLocalMaximumCWT(wCoefs)
#' ridgeList <- getRidge(localMax)
#'
#' SNR.Th <- 3
#' majorPeakInfo <- identifyMajorPeaks(exampleMS, ridgeList, wCoefs, SNR.Th = SNR.Th)
#' ## Plot the identified peaks
#' peakIndex <- majorPeakInfo$peakIndex
#' plotPeak(exampleMS, peakIndex, main = paste("Identified peaks with SNR >", SNR.Th))
#'
identifyMajorPeaks <- function(ms, ridgeList, wCoefs, scales = as.numeric(colnames(wCoefs)), SNR.Th = 3, peakScaleRange = 5,
    ridgeLength = 32, nearbyPeak = FALSE, nearbyWinSize = ifelse(nearbyPeak, 150, 100), winSize.noise = 500, SNR.method = "quantile", minNoiseLevel = 0.001, excludeBoundariesSize = nearbyWinSize / 2) {
    if (is.null(scales)) {
        scales <- 1:ncol(wCoefs)
        colnames(wCoefs) <- scales
    } else if (is.character(scales)) {
        scales <- as.numeric(scales)
    }
    if (ridgeLength > max(scales)) ridgeLength <- max(scales)

    if (length(peakScaleRange) == 1) {
        peakScaleRange <- scales[scales >= peakScaleRange]
    } else {
        peakScaleRange <- scales[scales >= peakScaleRange[1] & scales <= peakScaleRange[2]]
    }

    ## Limit the minNoiseLevel to avoid the case of very low noise level, e.g., smoothed spectrum
    if (minNoiseLevel >= 1) names(minNoiseLevel) <- "fixed"
    if (is.null(minNoiseLevel)) {
        minNoiseLevel <- 0
    } else { # By default the threshold is the ratio of the maximum coefficient
        if (is.null(names(minNoiseLevel))) {
            minNoiseLevel <- max(wCoefs) * minNoiseLevel
        } else if (names(minNoiseLevel) != "fixed") {
            minNoiseLevel <- max(wCoefs) * minNoiseLevel
        }
    }

    ## Get the peak values
    # mzInd <- as.numeric(names(ridgeList))
    ridgeLen <- sapply(ridgeList, length)
    ridgeName <- names(ridgeList)
    ridgeInfo <- matrix(as.numeric(unlist(strsplit(ridgeName, "_"))), nrow = 2)
    ridgeLevel <- ridgeInfo[1, ]
    # mzInd <- sapply(ridgeList, function(x) x[1])
    notnull <- sapply(ridgeList, function(x) {
        !is.null(x[1])
    }) # fixed by Steffen Neumann
    mzInd <- sapply(ridgeList[notnull], function(x) {
        x[1]
    }) # fixed by Steffen Neumann
    # mzInd <- ridgeInfo[2,]

    ## Reorder them by m/z index
    ord <- order(mzInd)
    ridgeName <- ridgeName[ord]
    ridgeLen <- ridgeLen[ord]
    ridgeLevel <- ridgeLevel[ord]
    ridgeList <- ridgeList[ord]
    mzInd <- mzInd[ord]

    peakScale <- NULL
    peakCenterInd <- NULL
    peakValue <- NULL
    # ridgeValue <- NULL
    ## Get the ridge values within the provided peakScaleRange
    for (i in 1:length(ridgeList)) {
        ridge.i <- ridgeList[[i]]
        level.i <- ridgeLevel[i]
        levels.i <- level.i:(level.i + ridgeLen[i] - 1)
        scales.i <- scales[levels.i]
        # Only keep the scales within the peakScaleRange
        selInd.i <- which(scales.i %in% peakScaleRange)
        if (length(selInd.i) == 0) {
            peakScale <- c(peakScale, scales.i[1])
            peakCenterInd <- c(peakCenterInd, ridge.i[1])
            peakValue <- c(peakValue, 0)
            next
        }

        levels.i <- levels.i[selInd.i]
        scales.i <- scales.i[selInd.i]
        ridge.i <- ridge.i[selInd.i]
        if (scales.i[1] == 0) {
            ind.i <- cbind(ridge.i[-1], levels.i[-1])
        } else {
            ind.i <- cbind(ridge.i, levels.i)
        }
        ridgeValue.i <- wCoefs[ind.i]
        maxInd.i <- which.max(ridgeValue.i)
        peakScale <- c(peakScale, scales.i[maxInd.i])
        peakCenterInd <- c(peakCenterInd, ridge.i[maxInd.i])
        peakValue <- c(peakValue, ridgeValue.i[maxInd.i])
        # ridgeValue <- c(ridgeValue, list(ridgeValue))
    }
    # names(ridgeValue) <- names(ridgeList)
    # ridgeLen <- get.ridgeLength(ridgeValue, Th=0.5)

    ## Compute SNR of each peak
    noise <- abs(wCoefs[, "1"])
    peakSNR <- NULL
    nMz <- nrow(wCoefs) # The length of ms signal

    for (k in 1:length(ridgeList)) {
        ind.k <- mzInd[k]
        start.k <- ifelse(ind.k - winSize.noise < 1, 1, ind.k - winSize.noise)
        end.k <- ifelse(ind.k + winSize.noise > nMz, nMz, ind.k + winSize.noise)
        ms.int <- ms[start.k:end.k] ## m/z intensity values in ind.k + /- winSize.noise (Added by Steffen Neumann)
        noiseLevel.k <- switch(SNR.method,
            quantile = stats::quantile(noise[start.k:end.k], probs = 0.95),
            sd = stats::sd(noise[start.k:end.k]),
            mad = stats::mad(noise[start.k:end.k], center = 0),
            data.mean = mean(ms.int), # (data.mean and data.mean.quant were added by Steffen Neumann)
            data.mean.quant = mean(ms.int[ms.int < stats::quantile(ms.int, probs = .95)]),
            stop("Invalid SNR.method. Please use one of 'quantile', 'sd', 'mad', 'data.mean', 'data.mean.quant'.")
        )
        ## Limit the minNoiseLevel to avoid the case of very low noise level, e.g., smoothed spectrum
        if (noiseLevel.k < minNoiseLevel) noiseLevel.k <- minNoiseLevel
        peakSNR <- c(peakSNR, peakValue[k] / noiseLevel.k)
    }

    ## Rule 1: ridge length should larger than a certain threshold
    # selInd1 <- (scales[ridgeLen] >= ridgeLength)
    selInd1 <- (scales[ridgeLevel + ridgeLen - 1] >= ridgeLength)

    ## In the case of nearbyPeak mode, it will include the nearby peaks within a certain range
    if (nearbyPeak) {
        selInd1 <- which(selInd1)
        index <- 1:length(mzInd)
        tempInd <- NULL
        for (ind.i in selInd1) {
            tempInd <- c(tempInd, index[mzInd >= mzInd[ind.i] - nearbyWinSize & mzInd <= mzInd[ind.i] + nearbyWinSize])
        }
        selInd1 <- (index %in% tempInd)
    }

    ## Rule 2: Based on the peak SNR
    if (SNR.Th > 0) {
    	selInd2 <- (peakSNR > SNR.Th)
    } else {
    	selInd2 <- rep(TRUE, length(peakSNR))
    }

    ## Because of the boundary effects,
    ## remove the peaks (half of the excludeBoundariesSize) at both ends of the signal profile if exists
    if (excludeBoundariesSize > 0) {
    	selInd3 <- !(mzInd %in% c(1:excludeBoundariesSize, (nrow(wCoefs) - excludeBoundariesSize + 1):nrow(wCoefs)))
    } else {
    	selInd3 <- rep(TRUE, length(mzInd))
    }

    ## combine SNR and peak length rule and other rules
    selInd <- (selInd1 & selInd2 & selInd3)

    names(peakSNR) <- names(peakScale) <- names(peakCenterInd) <- names(peakValue) <- names(mzInd) <- ridgeName

    list(
        peakIndex = mzInd[selInd],
        peakValue = peakValue,
        peakCenterIndex = peakCenterInd,
        peakSNR = peakSNR,
        peakScale = peakScale,
        potentialPeakIndex = mzInd[selInd1 & selInd3],
        allPeakIndex = mzInd
    )
}
