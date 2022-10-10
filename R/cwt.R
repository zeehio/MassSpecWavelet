#' Continuous Wavelet Transform (CWT)
#'
#' CWT(Continuous Wavelet Transform) with Mexican Hat wavelet (by default) to
#' match the peaks in Mass Spectrometry spectrum
#'
#' The default mother wavelet is a Mexican Hat evaluated in the \eqn{[-8,8]} range
#' using 1024 points (a step of roughly 1/64):
#'
#' \deqn{ \psi(x) = \frac{2}{\sqrt{3}} \pi^{-0.25} ( 1 - x^2 ) \exp{-x^2/2} }
#'
#' The \eqn{\sigma} of the mother Mexican Hat is of one x unit.
#'
#' The scaled wavelet is a downsampled version of the mother wavelet. The
#' `scale` determines how many samples per \eqn{x} unit are taken. For
#' instance, with the default Mexican Hat wavelet, a `scales = 2` will
#' evaluate the \eqn{[-8, 8]} range sampling twice per \eqn{x} unit, this is with a
#' step of 0.5. This generates a 33 points long scaled wavelet. Choosing this
#' type of scaling is convenient because the scaled wavelet becomes a wavelet
#' of \eqn{\sigma = `scales`} points. Using the default wavelet, a
#' `scales` smaller than 1 will show sampling issues, while a
#' `scales` larger than 64 will resample points from the original mother
#' wavelet. If you need to use `scales` larger than 64, consider providing
#' your own mother wavelet. See the examples.
#'
#' According to \doi{10.1063/1.3505103}, if
#' your spectrum has a gaussian peak shape of variance \eqn{m^2} points then
#' the scales range should cover \eqn{[1, 1.9 m]}. If your spectrum has a
#' Lorentzian peak shape of half-width-half-maximum \eqn{L} points then the
#' scales range should cover \eqn{[1, 2.9 L]}. They also suggest using a
#' \eqn{\log_{1.18}} spacing. Take these values as suggestions for your
#' analysis.
#'
#' @param ms Mass Spectrometry spectrum (a vector of MS intensities)
#' @param scales a vector represents the scales at which to perform CWT. See
#' the Details section. Additionally, a `prepared_wavelets` object
#' is also accepted (see [prepareWavelets()]).
#' @param wavelet The wavelet base, Mexican Hat by default. User can provide
#' wavelet `Psi(x)` as a form of two row matrix. The first row is the `x` value,
#' and the second row is `Psi(x)` corresponding to `x`.
#' @return The return is the 2-D CWT coefficient matrix, with column names as
#' the scale. Each column is the CWT coefficients at that scale. If the scales are
#' too big for the given signal, the returned matrix may include less columns than
#' the given scales.
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
    wavelets <- prepareWavelets(
        mslength = length(ms),
        scales = scales,
        wavelet = wavelet,
        wavelet_xlimit = 8,
        wavelet_length = 1024L,
        extendLengthScales = FALSE
    )
    scales <- wavelets$scales

    oldLen <- length(ms)
    
    # IF extendLengthScales is TRUE:
    # The new length is determined by the scales argument. See 
    # https://github.com/sneumann/xcms/issues/445 for more information.
    # extendLengthScales is now `FALSE` to preserve backwards compatibility,
    # but it should become `TRUE` in a future edition/version of the package
    if (wavelets$extendLengthScales) {
        ms <- extendLength(x = ms, addLength = (wavelets$padded_length - length(ms)), 
                               method = "open")
    } else {
        ms <- extendNBase(ms, nLevel = NULL, base = 2)
    }
    ms_fft <- stats::fft(ms)
    len <- length(ms)
    wCoefs <- matrix(NA_real_, nrow = oldLen, ncol = length(scales))

    for (i in seq_along(scales)) {
        scale.i <- scales[i]
        # Convolution:
        wCoefs.i <- 1 / sqrt(scale.i) * Re(stats::fft(ms_fft * wavelets$coefs[[i]], inverse = TRUE))/length(ms_fft)
        lenWave <- wavelets$lenWaves[i]
        ## Shift the position with half wavelet width
        wCoefs.i <- c(wCoefs.i[(len - floor(lenWave / 2) + 1):len], wCoefs.i[1:(len - floor(lenWave / 2))])
        wCoefs[, i] <- wCoefs.i[1:oldLen]
    }
    wCoefs <- wCoefs[,seq_len(i)]
    colnames(wCoefs) <- scales[seq_len(i)]
    return(wCoefs)
}

#' Prepare daughter wavelets for faster CWT
#' 
#' @param mslength Length of the signal to transform
#' @inheritParams cwt
#' @param wavelet_xlimit The mother wavelet will be evaluated between these limits. Ignored if `wavelet` is a matrix.
#' @param wavelet_length The number of points of the mother wavelet. Ignored if `wavelet` is a matrix
#' @param extendLengthScales A logical value. If the signal is too short, we may
#' need to pad it to convolve it with larger daughter wavelets. Set this to `TRUE` to let
#' scales be used to determine the padding length. It's set to `FALSE` by default
#' to preserve backwards compatibility.
#' @return A `prepared_wavelets` object.
#' @export
#' @examples 
#' x <- runif(2000)
#' scales <- c(1, 2, 4, 8)
#' prep_wavelets <- prepareWavelets(length(x), scales = scales)
#' wCoefs <- cwt(x, prep_wavelets)
#' @seealso cwt
prepareWavelets <- function(mslength, scales = c(1, seq(2, 30, 2), seq(32, 64, 4)), 
                             wavelet = "mexh", wavelet_xlimit = 8,
                             wavelet_length = 1024L, extendLengthScales = TRUE) {
    ## Check for the wavelet format
    if (inherits(scales, "prepared_wavelets")) {
        prepwavelet <- scales
        if (prepwavelet$mslength != mslength) {
            stop("The wavelets were not prepared for this ms length")
        }
        return(prepwavelet)
    }
    if (identical(wavelet, "mexh")) {
        psi_xval <- seq(-wavelet_xlimit, wavelet_xlimit, length.out = wavelet_length)
        psi <- mexh(psi_xval)
        # plot(psi_xval, psi)
    } else if (is.function(wavelet)) {
        psi_xval <- seq(-wavelet_xlimit, wavelet_xlimit, length.out = wavelet_length)
        psi <- wavelet(psi_xval)
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
    if (extendLengthScales) {
        len <- max(
            2^(ceiling(log2(max(scales)*(2*wavelet_xlimit)))),
            mslength
        )
    } else {
        len <- max(
            stats::nextn(mslength, 2),
            mslength
        )
    }
    
    psi_xval <- psi_xval - psi_xval[1]
    dxval <- psi_xval[2]
    xmax <- psi_xval[length(psi_xval)]
    prepared_wavelets <- vector("list", length = length(scales))
    lenWaves <- integer(length(scales))
    for (i in seq_along(scales)) {
        scale.i <- scales[i]
        f <- rep(0, len)
        j <- 1 + floor((0:(scale.i * xmax)) / (scale.i * dxval))
        if (length(j) == 1)
            j <- c(1, 1)
        lenWave <- length(j)
        f[1:lenWave] <- rev(psi[j]) - mean(psi[j])
        if (length(f) > len) {
            break
        }
        lenWaves[i] <- lenWave
        prepared.i <- Conj(stats::fft(f))
        prepared_wavelets[[i]] <- prepared.i
    }
    good_scales <- vapply(prepared_wavelets, function(x) !is.null(x), logical(1L))
    scales <- scales[good_scales]
    prepared_wavelets <- prepared_wavelets[good_scales]
    prepwavelet <- list(
        coefs = prepared_wavelets,
        scales = scales,
        mslength = mslength,
        padded_length = len,
        lenWaves = lenWaves,
        wavelet = wavelet,
        wavelet_xlimit = wavelet_xlimit,
        wavelet_length = wavelet_length,
        extendLengthScales = extendLengthScales
    )
    class(prepwavelet) <- "prepared_wavelets"
    prepwavelet
}


cwt_classic <- function(ms, scales=1, wavelet='mexh') {
        ## Check for the wavelet format
        if (identical(wavelet, 'mexh')) {
            psi_xval <- seq(-8, 8, length.out=1024)
            psi <- (2/sqrt(3) * pi^(-0.25)) * (1 - psi_xval^2) *exp(-psi_xval^2/2)
            #plot(psi_xval, psi)
        } else if (is.matrix(wavelet)) {
            if (nrow(wavelet) == 2) {
                psi_xval <- wavelet[1,]
                psi <- wavelet[2,]
            } else if (ncol(wavelet) == 2) {
                psi_xval <- wavelet[,1]
                psi <- wavelet[,2]
            } else {
                stop('Unsupported wavelet format!')
            }
        } else {
            stop('Unsupported wavelet!')
        }
        
        oldLen <- length(ms)
        ## To increase the computation effeciency of FFT, extend it as the power of 2
        ## because of a strange signal length 21577 makes the FFT very slow!
        ms <- extendNBase(ms, nLevel=NULL, base=2)
        len <- length(ms)
        nbscales <- length(scales)
        wCoefs <- NULL
        
        psi_xval <- psi_xval - psi_xval[1]
        dxval <- psi_xval[2]
        xmax  <- psi_xval[length(psi_xval)]
        for (i in 1:length(scales)) {
            scale.i <- scales[i]
            f <- rep(0, len)
            j <- 1 + floor((0:(scale.i * xmax))/(scale.i * dxval))
            if (length(j) == 1)		j <- c(1, 1)
            lenWave <- length(j)
            f[1:lenWave] <- rev(psi[j]) - mean(psi[j])
            if (length(f) > len) stop(paste('scale', scale.i, 'is too large!'))
            wCoefs.i <- 1/sqrt(scale.i) * convolve(ms, f)
            ## Shift the position with half wavelet width
            wCoefs.i <- c(wCoefs.i[(len-floor(lenWave/2) + 1) : len], wCoefs.i[1:(len-floor(lenWave/2))])
            wCoefs <- cbind(wCoefs, wCoefs.i)
        }
        if (length(scales) == 1) wCoefs <- matrix(wCoefs, ncol=1)
        colnames(wCoefs) <- scales
        wCoefs <- wCoefs[1:oldLen,,drop=FALSE]
        return(wCoefs)
    }


