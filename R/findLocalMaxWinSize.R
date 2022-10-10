#' Find local maxima and return the size of the window where they are maximum.
#' 
#' Compared to the rest of the package, this is a rather experimental function.
#' If you plan to use it or are interested in it, please open an issue at
#' https://github.com/zeehio/MassSpecWavelet/issues to show your interest.
#' 
#' @param x A numeric vector.
#' @param capWinSize the maximum window size to report. `NA` means unlimited.
#' @return An integer vector `y` of the same length as `x`. `y[i]` will be the
#' size of the largest window on `x` containing `x[i]` where:
#'  - `x[i]` is a local maximum or a center of a plateau
#'  - `x[i]` is not at a window border
#'  Optionally, if `capWinSize` is a positive integer, the maximum window size
#'  is capped to that value, to increase performance. Use this in case you just
#'  want to check if there exists a window of that size.
#'  @export
#'  @examples 
#'  x <- c(1, 2, 3, 2, 1)
#'  findLocalMaxWinSize(x)
findLocalMaxWinSize <- function(x, capWinSize = NA) {
    .Call("c_findLocalMaxWinSize", as.double(x), capWinSize = as.integer(capWinSize))
}