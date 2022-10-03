#' Plot the identified peaks over the spectrum
#'
#' Plot the identified peaks over the spectrum. The identified peaks are
#' returned by [peakDetectionCWT()] or
#' [identifyMajorPeaks()]
#'
#'
#' @param ms the MS spectrum
#' @param peakIndex m/z indexes of the identified peaks
#' @param mz m/z value correspond to m/z index
#' @param range the plot range of m/z value
#' @param method plot method of the identified peaks. method 'p' plot circles
#' on the peaks; method 'l' add vertical lines over the peaks.
#' @param main parameter of [plot()]
#' @param log parameter of [plot()]
#' @param \dots other parameters of [points()]
#' @return No value is returned; this function is called for its side effects
#' (plot).
#' @author Pan Du
#' @export
#' @seealso [peakDetectionCWT()], [identifyMajorPeaks()]
#' @keywords hplot
#' @examples
#'
#' data(exampleMS)
#' SNR.Th <- 3
#' peakInfo <- peakDetectionCWT(exampleMS, SNR.Th = SNR.Th)
#' majorPeakInfo <- peakInfo$majorPeakInfo
#' peakIndex <- majorPeakInfo$peakIndex
#' plotPeak(exampleMS, peakIndex, main = paste("Identified peaks with SNR >", SNR.Th))
#'
plotPeak <- function(ms, peakIndex = NULL, mz = 1:length(ms), range = c(min(mz), max(mz)), method = c("p", "l"), main = NULL, log = "", ...) {

    ## Check range parameter
    range <- sort(range)
    if (length(range) == 1) range <- c(range, max(mz))
    if (range[1] < min(mz)) range[1] <- min(mz)
    if (range[2] > max(mz)) range[2] <- max(mz)

    xlab <- ifelse(max(mz) == length(ms), "m/z index", "m/z value")

    selInd <- (mz > range[1] & mz < range[2])
    plot(mz[selInd], ms[selInd], type = "l", xlab = xlab, ylab = "Intensity", main = main, log = log)
    if (!is.null(peakIndex)) {
        if (method[1] == "p") {
            graphics::points(mz[peakIndex], ms[peakIndex], col = "red", ...)
        } else {
            graphics::abline(v = mz[peakIndex], col = "blue", ...)
        }
    }
}
