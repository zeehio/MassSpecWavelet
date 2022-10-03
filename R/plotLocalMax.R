#' Plot the local maximum matrix
#'
#' Plot the local maximum matrix of 2-D CWT coefficients returned by
#' \code{\link{getLocalMaximumCWT}}
#'
#'
#' @param localMax local maximum matrix of 2-D CWT coefficients returned by
#' \code{\link{getLocalMaximumCWT}}
#' @param wCoefs 2-D CWT coefficients
#' @param range plot range of m/z index
#' @param colorMap the colormap used in plotting the points
#' @param main parameter of \code{\link{plot}}
#' @param cex parameter of \code{\link{plot}}
#' @param pch parameter of \code{\link{plot}}
#' @param \dots other parameters of \code{\link{points}}
#' @return No value is returned; this function is called for its side effects
#' (plot).
#' @author Pan Du
#' @seealso \code{\link{getLocalMaximumCWT}}
#' @keywords hplot
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
plotLocalMax <- function(localMax, wCoefs = NULL, range = c(1, nrow(localMax)), colorMap = "RYB", main = NULL, cex = 3, pch = ".", ...) {
    if (colorMap == "RYB") {
        rgb.palette <- grDevices::colorRampPalette(c("red", "orange", "blue"), space = "rgb")
        colorMap <- rgb.palette(255)
        # colorMap <- heat.colors(255)
    }

    ## Check range parameter
    range <- sort(range, decreasing = FALSE)
    if (length(range) == 1) range <- c(range, ncol(localMax))
    if (range[1] < 1) range[1] <- 1
    if (range[2] > nrow(localMax)) range[2] <- ncol(localMax)

    test <- localMax[range[1]:range[2], ]
    tR <- range[1] - 1 + row(test)
    tC <- col(test)
    ind <- which(test > 0)
    if (is.null(wCoefs)) {
        col <- "blue"
    } else {
        selRidgeValue <- wCoefs[range[1]:range[2], ][ind]
        maxSelV <- log(1 + max(selRidgeValue))
        col <- colorMap[1 + log(1 + abs(selRidgeValue)) * 255 / maxSelV]
    }

    plot(seq(range[1], range[2], length = ncol(localMax)), 1:ncol(localMax), type = "n", xlab = "m/z index", ylab = "CWT coefficient scale", main = main)
    ind <- which(localMax[range[1]:range[2], ] > 0)
    graphics::points(tR[ind], tC[ind], col = col, cex = cex, pch = pch, ...)
}
