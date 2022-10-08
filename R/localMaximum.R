#' Identify local maximum within a slide window.
#'
#' Find local maximum by transform the vector as matrix, then get the the
#' maximum of each column. This operation is performed twice with vector
#' shifted half of the `winSize`.
#'
#' Instead of find the local maximum by a slide window, which slide all
#' possible positions, we find local maximum by transform the vector as matrix,
#' then get the the maximum of each column. This operation is performed twice
#' with vector shifted half of the `winSize`. The main purpose of this is to
#' increase the efficiency of the algorithm.
#'
#' @param x a vector represents a signal profile
#' @param winSize the slide window size, 5 by default.
#' @return Return a vector with the same length of the input x. The position of
#' local maximum is set as 1L, 0L else where.
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
    algo <- getOption("MassSpecWavelet.localMaximum.algorithm", "new")
    if (!algo %in% c("new", "classic")) {
        warning('Invalid algorithm "', algo, '". Use either "new" or "classic". Assuming "classic".')
        algo <- "classic"
    }
    if (algo == "new") {
        local_max <- findLocalMaxWinSize(x, capWinSize = winSize)
        return(as.integer(local_max >= winSize))
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
