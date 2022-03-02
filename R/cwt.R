
get_scaled_wavelet <- local({
	# We use a cache of computed wavelets to avoid recomputations
	cache <- list()
	function(wavelet, scale, target_len) {
		cache_key <- NULL
		if (identical(wavelet, 'clear_cache')) {
			cache[] <- NULL
			return(NULL)
        }
		## Check for the wavelet format
		if (identical(wavelet, 'mexh')) {
			cache_key <- paste0(wavelet, "_", as.character(scale), "_", as.character(target_len))
			if (cache_key %in% names(cache)) {
				return(cache[[cache_key]])
			}
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

		psi_xval <- psi_xval - psi_xval[1]
		dxval <- psi_xval[2]
		xmax  <- psi_xval[length(psi_xval)]

		f <- rep(0, target_len)
        j <- 1 + floor((0:(scale * xmax))/(scale * dxval))
        if (length(j) == 1)		j <- c(1, 1)
		lenWave <- length(j)
        f[1:lenWave] <- rev(psi[j]) - mean(psi[j])

		conj_fft_f <- Conj(fft(f))
		out <- list(
			conj_fft_f = conj_fft_f,
			lenWave = lenWave,
			is_wavelet_real = is.numeric(psi)
		)
		if (!is.null(cache_key)) {
			cache[[cache_key]] <<- out
		}
		out
	}
})

"cwt" <-
function(ms, scales=1, wavelet='mexh') {
	if (identical(wavelet, 'clear_cache')) {
		return(get_scaled_wavelet(wavelet, 0, 0))
	}
	oldLen <- length(ms)
	## To increase the computation effeciency of FFT, extend it as the power of 2
	## because of a strange signal length 21577 makes the FFT very slow!
	ms <- extendNBase(ms, nLevel=NULL, base=2)
	ms_fft <- fft(ms)
	len <- length(ms)
    nbscales <- length(scales)
    wCoefs <- NULL


    for (i in 1:length(scales)) {
		scale.i <- scales[i]
		wv <- get_scaled_wavelet(wavelet, scale.i, len)
		if (length(wv$conj_fft_f) > len) stop(paste('scale', scale.i, 'is too large!'))
		convolved <- fft(ms_fft * wv$conj_fft_f, inverse = TRUE)
		if (is.numeric(ms) & wv$is_wavelet_real) {
			convolved <- Re(convolved)
		}
		wCoefs.i <- 1/sqrt(scale.i) * convolved/len
		## Shift the position with half wavelet width
		wCoefs.i <- c(wCoefs.i[(len-floor(wv$lenWave/2) + 1) : len], wCoefs.i[1:(len-floor(wv$lenWave/2))])
		wCoefs <- cbind(wCoefs, wCoefs.i)
    }
	if (length(scales) == 1) wCoefs <- matrix(wCoefs, ncol=1)
	colnames(wCoefs) <- scales
	wCoefs <- wCoefs[1:oldLen,,drop=FALSE]
	return(wCoefs)
}

