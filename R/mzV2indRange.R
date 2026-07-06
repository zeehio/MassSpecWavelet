#' Match m/z value to m/z index with a certain error range
#'
#' Match m/z value to m/z index with a certain error range
#'
#' **Deprecated**: this function appears to be unused, both internally in
#' MassSpecWavelet and by known downstream packages, and is scheduled for
#' removal in a future release. If you rely on it, please open an issue at
#' <https://github.com/zeehio/MassSpecWavelet/issues> within the next year.
#'
#' @param mzV a vector of m/z value
#' @param error error range
#' @return return a vector of sorted m/z indexes
#' @author Pan Du
#' @seealso [mzInd2vRange()]
#' @keywords methods
mzV2indRange <- function(mzV, error = 0.003) {
    .Deprecated(msg = paste(
        "mzV2indRange() appears to be unused, both internally in MassSpecWavelet",
        "and by known downstream packages, and is scheduled for removal in a",
        "future release. If you rely on it, please open an issue at",
        "https://github.com/zeehio/MassSpecWavelet/issues within the next year."
    ))
    mzIndR <- NULL
    for (i in 1:length(mzV)) {
        # suppressWarnings(): avoid repeating u2i()'s own deprecation warning
        # once per loop iteration; the warning above already covers this call.
        from.i <- round(suppressWarnings(u2i(mzV[i] * (1 - error))))
        to.i <- round(suppressWarnings(u2i(mzV[i] * (1 + error))))
        mzIndR <- c(mzIndR, from.i:to.i)
    }
    mzIndR <- sort(unique(mzIndR))
    return(mzIndR)
}
