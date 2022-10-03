#' Get the CWT coefficient values corresponding to the peak ridge
#'
#' Get the CWT coefficient values corresponding to the peak ridge
#'
#'
#' @param ridgeList a list of ridge lines
#' @param wCoefs 2-D CWT coefficients
#' @param skip the CWT scale level to be skipped, by default the 0 scale level
#' (raw spectrum) is skipped.
#' @return A list of ridge values corresponding to the input ridgeList.
#' @author Pan Du
#' @keywords methods
getRidgeValue <- function(ridgeList, wCoefs, skip = 0) {
    ridgeLen <- vapply(ridgeList, length, integer(1L))
    ridgeName <- names(ridgeList)
    ridgeInfo <- matrix(as.numeric(unlist(strsplit(ridgeName, "_"))), nrow = 2)
    ridgeLevel <- ridgeInfo[1, ]
    mzInd <- ridgeInfo[2, ]

    scales <- attr(ridgeList, "scales")
    if (is.null(scales)) scales <- 0:max(ridgeLen)
    if (length(scales) != ncol(wCoefs)) stop('The length of "scales" does not match with "wCoefs"!')
    skip <- as.character(skip)

    ridgeValue <- NULL
    ridgeMax <- NULL
    for (i in 1:length(mzInd)) {
        ridge.i <- ridgeList[[i]]
        level.i <- ridgeLevel[i]
        levels.i <- level.i:(level.i + ridgeLen[i] - 1)
        scales.i <- scales[levels.i]
        removeInd.i <- which(scales.i == skip)
        if (length(removeInd.i) > 0) {
            levels.i <- levels.i[-removeInd.i]
            ridge.i <- ridge.i[-removeInd.i]
            scales.i <- scales.i[-removeInd.i]
        }
        coefInd.i <- cbind(ridge.i, levels.i)
        ridgeValue.i <- wCoefs[coefInd.i]
        names(ridgeValue.i) <- scales.i
        ridgeValue <- c(ridgeValue, list(ridgeValue.i))
        ridgeMax <- c(ridgeMax, max(ridgeValue.i))
    }
    names(ridgeValue) <- names(ridgeMax) <- ridgeName
    attr(ridgeValue, "ridgeMax") <- ridgeMax
    return(ridgeValue)
}
