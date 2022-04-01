"plotPeak" <-
function(ms, peakIndex=NULL, mz=1:length(ms), range=c(min(mz), max(mz)), method=c('p', 'l'), main=NULL, log='', ... ) {
	
	## Check range parameter
	range <- sort(range)
	if (length(range) == 1) range <- c(range, max(mz))
	if (range[1] < min(mz)) range[1] <- min(mz)
	if (range[2] > max(mz)) range[2] <- max(mz)
	
	xlab <- ifelse (max(mz) == length(ms), 'm/z index', 'm/z value')
	
	selInd <- (mz > range[1] & mz < range[2])
	plot(mz[selInd], ms[selInd], type='l', xlab=xlab, ylab='Intensity', main=main, log=log)
	if (!is.null(peakIndex)) {
		if (method[1] == 'p') {
			graphics::points(mz[peakIndex], ms[peakIndex], col='red', ...)
		} else {
			graphics::abline(v=mz[peakIndex], col='blue',...)
		}
	}
}

