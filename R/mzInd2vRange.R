#' Match m/z index to m/z value with a certain error range
#'
#' Match m/z index to m/z value with a certain error range
#'
#'
#' @param mzInd a vector of m/z index
#' @param error error range
#' @return return a vector of sorted m/z values
#' @author Pan Du
#' @seealso \code{\link{mzV2indRange}}
#' @keywords methods
mzInd2vRange <- function(mzInd, error = 0.003) {
    mzVR <- NULL
    for (i in 1:length(mzInd)) {
        from.i <- round(i2u(mzInd[i]) * (1 - error))
        to.i <- round(i2u(mzInd[i]) * (1 + error))
        mzVR <- c(mzVR, from.i:to.i)
    }
    mzVR <- sort(unique(mzVR))
    return(mzVR)
}
