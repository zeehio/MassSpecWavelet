#' Identify the local maximum of each column in 2-D CWT coefficients matrix
#'
#' Identify the local maximum of each column in 2-D CWT coefficients matrix by
#' using a slide window. The size of slide window linearly changes from the
#' coarse scale (bigger window size) to detail scale. The scale of CWT
#' increases with the column index.
#'
#'
#' @param wCoefs 2-D CWT coefficients, each column corresponding to CWT
#' coefficient at one scale. The column name is the scale.
#' @param minWinSize The minimum slide window size used.
#' @param amp.Th The minimum peak amplitude.
#' @param is_amp_thres_relative Whether `amp.Th` is given relative to
#' `max(wCoefs)`.
#' @return return a matrix with same dimension as CWT coefficient matrix,
#' wCoefs. The local maxima are marked as 1, others are 0.
#' @author Pan Du
#' @seealso [localMaximum()]
#' @keywords methods
#' @export
#' @examples
#'
#' data(exampleMS)
#' scales <- seq(1, 64, 3)
#' wCoefs <- cwt(exampleMS[5000:11000], scales = scales, wavelet = "mexh")
#'
#' localMax <- getLocalMaximumCWT(wCoefs)
#' plotLocalMax(localMax)
#'
getLocalMaximumCWT <- function(wCoefs, minWinSize = 5, amp.Th = 0, is_amp_thres_relative = FALSE) {
    localMax <- NULL
    scales <- as.numeric(colnames(wCoefs))

    if (isTRUE(is_amp_thres_relative)) {
        amp.Th <- max(wCoefs) * amp.Th
    }

    for (i in seq_along(scales)) {
        scale.i <- scales[i]
        winSize.i <- scale.i * 2 + 1
        if (winSize.i < minWinSize) {
            winSize.i <- minWinSize
        }
        temp <- localMaximum(wCoefs[, i], winSize.i)
        localMax <- cbind(localMax, temp)
    }
    # Set the values less than peak threshold as 0
    localMax[wCoefs < amp.Th] <- 0L
    colnames(localMax) <- colnames(wCoefs)
    rownames(localMax) <- rownames(wCoefs)
    return(localMax)
}
