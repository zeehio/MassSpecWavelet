#' Match m/z index to m/z value with a certain error range
#'
#' Match m/z index to m/z value with a certain error range
#'
#' **Deprecated**: this function appears to be unused, both internally in
#' MassSpecWavelet and by known downstream packages, and is scheduled for
#' removal in a future release. If you rely on it, please open an issue at
#' <https://github.com/zeehio/MassSpecWavelet/issues> within the next year.
#'
#' @param mzInd a vector of m/z index
#' @param error error range
#' @return return a vector of sorted m/z values
#' @author Pan Du
#' @seealso [mzV2indRange()]
#' @keywords methods
mzInd2vRange <- function(mzInd, error = 0.003) {
    .Deprecated(msg = paste(
        "mzInd2vRange() appears to be unused, both internally in MassSpecWavelet",
        "and by known downstream packages, and is scheduled for removal in a",
        "future release. If you rely on it, please open an issue at",
        "https://github.com/zeehio/MassSpecWavelet/issues within the next year."
    ))
    mzVR <- NULL
    for (i in 1:length(mzInd)) {
        # suppressWarnings(): avoid repeating i2u()'s own deprecation warning
        # once per loop iteration; the warning above already covers this call.
        from.i <- round(suppressWarnings(i2u(mzInd[i])) * (1 - error))
        to.i <- round(suppressWarnings(i2u(mzInd[i])) * (1 + error))
        mzVR <- c(mzVR, from.i:to.i)
    }
    mzVR <- sort(unique(mzVR))
    return(mzVR)
}
