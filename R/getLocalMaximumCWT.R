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
#' @param isAmpThreshRelative Whether `amp.Th` is given relative to
#' `max(wCoefs)`.
#' @param exclude0scaleAmpThresh When computing the relative `amp.Th`, if
#' this is set to `TRUE`, the `amp.Th` will exclude the zero-th scale from the
#' `max(wCoefs)`. The zero-th scale corresponds to the original signal, that may
#' have a much larger baseline than the wavelet coefficients and can distort the
#' threshold calculation. The default is `FALSE` to preserve backwards compatibility.
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
getLocalMaximumCWT <- function(wCoefs, minWinSize = 5, amp.Th = 0, isAmpThreshRelative = FALSE, exclude0scaleAmpThresh = FALSE) {
    localMax <- matrix(NA_integer_, nrow = nrow(wCoefs), ncol = ncol(wCoefs))
    scales <- as.numeric(colnames(wCoefs))
    
    if (isTRUE(isAmpThreshRelative)) {
        if (isTRUE(exclude0scaleAmpThresh) && isTRUE("0" %in% colnames(wCoefs))) {
            amp.Th <- max(wCoefs[,colnames(wCoefs) != "0", drop = FALSE]) * amp.Th
        } else {
            amp.Th <- max(wCoefs) * amp.Th
        }
    }
    
    for (i in seq_along(scales)) {
        scale.i <- scales[i]
        winSize.i <- scale.i * 2 + 1
        if (winSize.i < minWinSize) {
            winSize.i <- minWinSize
        }
        localMax[, i] <- localMaximum(wCoefs[, i], winSize.i)
    }
    # Set the values less than peak threshold as 0
    localMax[wCoefs < amp.Th] <- 0L
    colnames(localMax) <- colnames(wCoefs)
    rownames(localMax) <- rownames(wCoefs)
    return(localMax)
}
