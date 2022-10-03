#' Match m/z value to m/z index with a certain error range
#'
#' Match m/z value to m/z index with a certain error range
#'
#'
#' @param mzV a vector of m/z value
#' @param error error range
#' @return return a vector of sorted m/z indexes
#' @author Pan Du
#' @seealso [mzInd2vRange()]
#' @keywords methods
mzV2indRange <- function(mzV, error = 0.003) {
    mzIndR <- NULL
    for (i in 1:length(mzV)) {
        from.i <- round(u2i(mzV[i] * (1 - error)))
        to.i <- round(u2i(mzV[i] * (1 + error)))
        mzIndR <- c(mzIndR, from.i:to.i)
    }
    mzIndR <- sort(unique(mzIndR))
    return(mzIndR)
}
