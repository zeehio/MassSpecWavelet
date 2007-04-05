"getRidgeLength" <-
function(ridgeList, Th=0.5) {

	ridgeLen <- sapply(ridgeList, function(x) {
		m <- max(x)
		l <- length(x)
		while(l > 1) {
			if (x[l] > m * Th) break
			l <- l - 1
		}
		return(l)
	} )
	return(ridgeLen)
}

