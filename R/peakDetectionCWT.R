"peakDetectionCWT" <-
function(ms, scales=c(1, seq(2,30,2),seq(32, 64, 4)), SNR.Th=3, nearbyPeak=TRUE, peakScaleRange=5, 
		amp.Th=0.01, minNoiseLevel=amp.Th/SNR.Th, ridgeLength=24, tuneIn=FALSE, ... ) {

	if (minNoiseLevel > 1)  names(minNoiseLevel) <- 'fixed' 
	## Perform Continuous Wavelet Transform
	wCoefs <- cwt(ms, scales=scales, wavelet='mexh')

	## Attach the raw data as the zero level of decomposition
	wCoefs <- cbind(as.vector(ms), wCoefs)
	colnames(wCoefs) <- c(0, scales)

	##-----------------------------------------
	## Identify the local maximum by using a slide window
	## The size of slide window changes over different levels, with the coarse level have bigger window size
	if (is.null(amp.Th)) {
		amp.Th <- 0 
	} else {
		if (is.null(names(amp.Th))) {
			amp.Th <- max(wCoefs) * amp.Th 
		} else if (names(amp.Th) != 'fixed') {
			amp.Th <- max(wCoefs) * amp.Th
		}
	} 
	localMax <- getLocalMaximumCWT(wCoefs, amp.Th=amp.Th, ...)
	colnames(localMax) <- colnames(wCoefs)

	##-----------------------------------------
	## Indentify the ridges from coarse level to more detailed levels
	ridgeList <- getRidge(localMax, gapTh=3, skip=2, ...)

	##-----------------------------------------
	## Indentify the major peaks and their nearby peaks 
	majorPeakInfo <- identifyMajorPeaks(ridgeList, wCoefs, SNR.Th=SNR.Th, peakScaleRange=peakScaleRange, 
			nearbyPeak=nearbyPeak, minNoiseLevel=minNoiseLevel, ridgeLength=ridgeLength, ...)

	if (tuneIn) {
		refinedPeakInfo <- tuneInPeakInfo(ms, majorPeakInfo)
		return(list(majorPeakInfo=refinedPeakInfo, ridgeList=ridgeList, localMax=localMax, wCoefs=wCoefs[, -1], oldPeakInfo=majorPeakInfo))		
	} else {
		return(list(majorPeakInfo=majorPeakInfo, ridgeList=ridgeList, localMax=localMax, wCoefs=wCoefs[, -1]))		
	}
}

