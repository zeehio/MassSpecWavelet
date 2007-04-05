"mzInd2vRange" <-
function(mzInd, error=0.003) {

	mzVR <- NULL
	for (i in 1:length(mzInd)) {
		from.i <- round(i2u(mzInd[i]) * (1 - error))
		to.i <- round(i2u(mzInd[i]) * (1 + error))
		mzVR <- c(mzVR, from.i:to.i)
	}
	mzVR <- sort(unique(mzVR))
	return(mzVR)
}

