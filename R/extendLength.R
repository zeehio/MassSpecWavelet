#' Extend the length of a signal or matrix
#'
#' Extend the length of a signal or matrix by row
#'
#'
#' @param x a vector or matrix with column with each column as a signal
#' @param addLength the length to be extended
#' @param method three methods available, c("reflection", "open", "circular").
#' By default, it is "reflection".
#' @param direction three options available: c("right", "left", "both")
#' @return return the extended vector or matrix.
#' @author Pan Du
#' @seealso [extendNBase()]
#' @keywords methods
#' @examples
#'
#' # a = matrix(rnorm(9), 3)
#' # extendLength(a, 3, direction='right') 	## not exposed function
#'
extendLength <- function(x, addLength = NULL, method = c("reflection", "open", "circular"), direction = c("right", "left", "both")) {
    if (is.null(addLength)) stop("Please provide the length to be added!")
    if (!is.matrix(x)) x <- matrix(x, ncol = 1)
    method <- match.arg(method)
    direction <- match.arg(direction)

    nR <- nrow(x)
    nR1 <- nR + addLength
    if (direction == "both") {
        left <- right <- addLength
    } else if (direction == "right") {
        left <- 0
        right <- addLength
    } else if (direction == "left") {
        left <- addLength
        right <- 0
    }

    if (right > 0) {
        x <- switch(method,
            reflection = rbind(x, x[nR:(2 * nR - nR1 + 1), , drop = FALSE]),
            open = rbind(x, matrix(rep(x[nR, ], addLength), ncol = ncol(x), byrow = TRUE)),
            circular = rbind(x, x[1:(nR1 - nR), , drop = FALSE])
        )
    }

    if (left > 0) {
        x <- switch(method,
            reflection = rbind(x[addLength:1, , drop = FALSE], x),
            open = rbind(matrix(rep(x[1, ], addLength), ncol = ncol(x), byrow = TRUE), x),
            circular = rbind(x[(2 * nR - nR1 + 1):nR, , drop = FALSE], x)
        )
    }
    if (ncol(x) == 1) x <- as.vector(x)

    return(x)
}
