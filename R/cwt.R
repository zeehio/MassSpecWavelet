#' Continuous Wavelet Transform (CWT)
#'
#' CWT(Continuous Wavelet Transform) with Mexican Hat wavelet (by default) to
#' match the peaks in Mass Spectrometry spectrum
#'
#' The default mother wavelet is a Mexican Hat evaluated in the [-8,8] range
#' using 1024 points (a step of roughly 1/64):
#'
#' \deqn{ \psi(x) = \frac{2}{\sqrt{3}} \pi^{-0.25} ( 1 - x^2 ) \exp{-x^2/2} }
#'
#' The \eqn{\sigma} of the mother Mexican Hat is of one x unit.
#'
#' The scaled wavelet is a downsampled version of the mother wavelet. The
#' `scale` determines how many samples per \eqn{x} unit are taken. For
#' instance, with the default Mexican Hat wavelet, a `scales = 2` will
#' evaluate the [-8, 8] range sampling twice per \eqn{x} unit, this is with a
#' step of 0.5. This generates a 33 points long scaled wavelet. Choosing this
#' type of scaling is convenient because the scaled wavelet becomes a wavelet
#' of \eqn{\sigma = `scales`} points. Using the default wavelet, a
#' `scales` smaller than 1 will show sampling issues, while a
#' `scales` larger than 64 will resample points from the original mother
#' wavelet. If you need to use `scales` larger than 64, consider providing
#' your own mother wavelet. See the examples.
#'
#' According to c("\\Sexpr[results=rd]{tools:::Rd_expr_doi(\"#1\")}",
#' "10.1063/1.3505103")\Sexpr{tools:::Rd_expr_doi("10.1063/1.3505103")}, if
#' your spectrum has a gaussian peak shape of variance \eqn{m^2} points then
#' the scales range should cover \eqn{[1, 1.9 m]}. If your spectrum has a
#' Lorentzian peak shape of half-width-half-maximum \eqn{L} points then the
#' scales range should cover \eqn{[1, 2.9 L]}. They also suggest using a
#' \eqn{\log_{1.18}} spacing. Take these values as suggestions for your
#' analysis.
#'
#' @param ms Mass Spectrometry spectrum (a vector of MS intensities)
#' @param scales a vector represents the scales at which to perform CWT. See
#' the Details section
#' @param wavelet The wavelet base, Mexican Hat by default. User can provide
#' wavelet Psi(x) as a form of two row matrix. The first row is the x value,
#' and the second row is Psi(x) corresponding to x.
#' @return The return is the 2-D CWT coefficient matrix, with column names as
#' the scale. Each column is the CWT coefficients at that scale.
#' @author Pan Du, Simon Lin
#' @keywords methods
#' @examples
#'
#' data(exampleMS)
#' scales <- seq(1, 64, 3)
#' wCoefs <- cwt(exampleMS[5000:11000], scales = scales, wavelet = "mexh")
#'
#' ## Plot the 2-D CWT coefficients as image (It may take a while!)
#' xTickInterval <- 1000
#' image(5000:11000, scales, wCoefs,
#'     col = terrain.colors(256), axes = FALSE,
#'     xlab = "m/z index", ylab = "CWT coefficient scale", main = "CWT coefficients"
#' )
#' axis(1, at = seq(5000, 11000, by = xTickInterval))
#' axis(2, at = c(1, seq(10, 64, by = 10)))
#' box()
#'
#' ## Provide a larger wavelet:
#' scales <- c(seq(1, 64, 3), seq(72, 128, 8))
#' psi_xval <- seq(-8, 8, length.out = 4096)
#' psi <- (2 / sqrt(3) * pi^(-0.25)) * (1 - psi_xval^2) * exp(-psi_xval^2 / 2)
#' wCoefs <- cwt(exampleMS[5000:11000], scales = scales, wavelet = rbind(psi_xval, psi))
#' xTickInterval <- 1000
#' image(5000:11000, scales, wCoefs,
#'     col = terrain.colors(256), axes = FALSE,
#'     xlab = "m/z index", ylab = "CWT coefficient scale", main = "CWT coefficients"
#' )
#' axis(1, at = seq(5000, 11000, by = xTickInterval))
#' axis(2, at = c(1, seq(10, 128, by = 10)))
#' box()
#'
#' ## Custom log1.18 spaced scales:
#' get_scales <- function(scale_min, scale_max, num_scales) {
#'     (seq(0, 1, length.out = num_scales)^1.18) * (scale_max - scale_min) + scale_min
#' }
#' scales <- get_scales(scale_min = 1, scale_max = 64, num_scales = 32)
#' print(scales)
#' # For detecting a gaussian peak of 10 points:
#' gaussian_peak_sigma <- 10 # points
#' scales <- get_scales(scale_min = 1, scale_max = 1.9 * gaussian_peak_sigma, num_scales = 32)
#' print(scales)
#' # For detecting a lorentzian peak of 10 points:
#' lorentzian_peak_gamma <- 10 # points
#' scales <- get_scales(scale_min = 1, scale_max = 2.9 * lorentzian_peak_gamma, num_scales = 32)
#' print(scales)
#'
#' @export
cwt <- function(ms, scales = 1, wavelet = "mexh") {
    ## Check for the wavelet format
    if (identical(wavelet, "mexh")) {
        psi_xval <- seq(-8, 8, length.out = 1024)
        psi <- (2 / sqrt(3) * pi^(-0.25)) * (1 - psi_xval^2) * exp(-psi_xval^2 / 2)
        # plot(psi_xval, psi)
    } else if (is.matrix(wavelet)) {
        if (nrow(wavelet) == 2) {
            psi_xval <- wavelet[1, ]
            psi <- wavelet[2, ]
        } else if (ncol(wavelet) == 2) {
            psi_xval <- wavelet[, 1]
            psi <- wavelet[, 2]
        } else {
            stop("Unsupported wavelet format!")
        }
    } else {
        stop("Unsupported wavelet!")
    }

    oldLen <- length(ms)
    ## To increase the computation effeciency of FFT, extend it as the power of 2
    ## because of a strange signal length 21577 makes the FFT very slow!
    ms <- extendNBase(ms, nLevel = NULL, base = 2)
    len <- length(ms)
    nbscales <- length(scales)
    wCoefs <- NULL

    psi_xval <- psi_xval - psi_xval[1]
    dxval <- psi_xval[2]
    xmax <- psi_xval[length(psi_xval)]
    for (i in 1:length(scales)) {
        scale.i <- scales[i]
        f <- rep(0, len)
        j <- 1 + floor((0:(scale.i * xmax)) / (scale.i * dxval))
        if (length(j) == 1) j <- c(1, 1)
        lenWave <- length(j)
        f[1:lenWave] <- rev(psi[j]) - mean(psi[j])
        if (length(f) > len) stop("scale ", scale.i, " is too large!")
        wCoefs.i <- 1 / sqrt(scale.i) * stats::convolve(ms, f)
        ## Shift the position with half wavelet width
        wCoefs.i <- c(wCoefs.i[(len - floor(lenWave / 2) + 1):len], wCoefs.i[1:(len - floor(lenWave / 2))])
        wCoefs <- cbind(wCoefs, wCoefs.i)
    }
    if (length(scales) == 1) wCoefs <- matrix(wCoefs, ncol = 1)
    colnames(wCoefs) <- scales
    wCoefs <- wCoefs[1:oldLen, , drop = FALSE]
    return(wCoefs)
}
