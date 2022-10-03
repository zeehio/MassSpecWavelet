#' Estimate the length of the ridge
#'
#' Estimate the length of the ridge line, which is composed of local maxima at
#' adjacent CWT scales. The ridge line is cut off at the end point, whose
#' amplitude divided by the maximum ridge amplitude is larger than the cutoff
#' amplitude ratio threshold (0.5 by default).
#'
#'
#' @param ridgeList a list of identified ridges
#' @param Th the cutoff amplitude ratio (the amplitude divided by the maximum
#' amplitude of the ridge) threshold of the ridge line end.
#' @return a vector of estimated ridge length
#' @author Pan Du
#' @keywords methods
#' @export
#' @examples 
#' stopifnot(getRidgeLength(list(c(5,4,3,2,1), c(5,3,1))) == c(3,2))
getRidgeLength <- function(ridgeList, Th = 0.5) {
    ridgeLen <- vapply(ridgeList, function(x) {
        m <- max(x)
        l <- length(x)
        while (l > 1) {
            if (x[l] > m * Th) break
            l <- l - 1
        }
        return(l)
    }, integer(1L))
    return(ridgeLen)
}
