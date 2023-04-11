#' Identify local maximum within a slide window.
#'
#' The simplest local maximum detection using a sliding window searches for
#' maxima in a window of a given size, and slides that window across the signal,
#' shifting it one position at a time.
#' 
#' The default implementation found here shifts the window by half of its size
#' instead of by one position at a time. This makes the implementation faster,
#' at the expense of not being able to detect peaks that are too close to each other,
#' if they appear in some positions with respect to the windows.
#' 
#' Additionally, this implementation removes all instances of peaks found at
#' a distance less than the window size
#' 
#' Experimentally, we are exploring other algorithms for local maxima detection.
#' These algorithms can be chosen setting the `"MassSpecWavelet.localMaximum.algorithm"`
#' option. See the `"Finding local maxima"` vignette for further details.
#' 
#' 
#' @param x a vector represents a signal profile
#' @param winSize the slide window size, 5 by default.
#' @return Return a vector with the same length of the input x. The position of
#' local maximum is set as `1L`, `0L` else where.
#' @author Pan Du and Sergio Oller
#' @seealso [getLocalMaximumCWT()]
#' @keywords methods
#' @export
#' @examples
#'
#' x <- rnorm(200)
#' lmax <- localMaximum(x, 5)
#' maxInd <- which(lmax > 0)
#' plot(x, type = "l")
#' points(maxInd, x[maxInd], col = "red")
#'
localMaximum <- function(x, winSize = 5) {
    algo <- getOption("MassSpecWavelet.localMaximum.algorithm", "faster")
    if (!algo %in% c("new", "classic", "faster")) {
        warning('Invalid algorithm "', algo, '". Use either "new", "faster" or "classic". Assuming "faster".')
        algo <- "faster"
    }
    if (algo == "new") {
        local_max <- findLocalMaxWinSize(x, capWinSize = winSize)
        localMax <- as.integer(local_max >= winSize)
        return(localMax)
    } else if (algo == "faster") {
        localMax <- localMaximumSlidingWindow(x, winSize)
        ## Check whether there is some local maxima have in between distance less than winSize
        maxInd <- which(localMax > 0)
        selInd <- which(diff(maxInd) < winSize)
        if (length(selInd) > 0) {
            selMaxInd1 <- maxInd[selInd]
            selMaxInd2 <- maxInd[selInd + 1L]
            temp <- x[selMaxInd1] - x[selMaxInd2]
            localMax[selMaxInd1[temp <= 0]] <- 0L
            localMax[selMaxInd2[temp > 0]] <- 0L
        }
        return(localMax)
    }
    len <- length(x)
    rNum <- ceiling(len / winSize)

    ## Transform the vector as a matrix with column length equals winSize
    ## and find the maximum position at each row.
    y <- matrix(c(x, rep(x[len], rNum * winSize - len)), nrow = winSize)
    y.maxInd <- apply(y, 2, which.max)
    ## Only keep the maximum value larger than the boundary values
    selInd <- which(apply(y, 2, function(x) max(x) > x[1] & max(x) > x[winSize]))

    ## keep the result
    localMax <- rep(0L, len)
    localMax[(selInd - 1) * winSize + y.maxInd[selInd]] <- 1L

    ## Shift the vector with winSize/2 and do the same operation
    shift <- floor(winSize / 2)
    rNum <- ceiling((len + shift) / winSize)
    y <- matrix(c(rep(x[1], shift), x, rep(x[len], rNum * winSize - len - shift)), nrow = winSize)
    y.maxInd <- apply(y, 2, which.max)
    ## Only keep the maximum value larger than the boundary values
    selInd <- which(apply(y, 2, function(x) max(x) > x[1] & max(x) > x[winSize]))
    localMax[(selInd - 1) * winSize + y.maxInd[selInd] - shift] <- 1L

    ## Check whether there is some local maxima have in between distance less than winSize
    maxInd <- which(localMax > 0)
    selInd <- which(diff(maxInd) < winSize)
    if (length(selInd) > 0) {
        selMaxInd1 <- maxInd[selInd]
        selMaxInd2 <- maxInd[selInd + 1L]
        temp <- x[selMaxInd1] - x[selMaxInd2]
        localMax[selMaxInd1[temp <= 0]] <- 0L
        localMax[selMaxInd2[temp > 0]] <- 0L
    }

    return(localMax)
}

localMaximumSlidingWindow <- function(x, winSize = 5L) {
    .Call(c_localMaximumSlidingWindow, as.double(x), as.integer(winSize))
}