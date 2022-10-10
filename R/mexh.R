#' The mexican hat function
#' 
#' \deqn{ \psi(x) = \frac{2}{\sqrt{3}} \pi^{-0.25} ( 1 - x^2 ) \exp{-x^2/2} }
#' 
#' @param x where to evaluate the mexican hat
#' @return A vector of the same length as `x` with the corresponding values
#' @examples
#' x <- seq(-8, 8, length.out = 256)
#' mexh(x)
#' @export
mexh <- function(x) {
    (2 / sqrt(3) * pi^(-0.25)) * (1 - x^2) * exp(-x^2 / 2)
}
