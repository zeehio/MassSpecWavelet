"smoothDWT" <-
function(ms, nLevel=6, wf="la8", localNoiseTh=seq(1, 0, by=-0.2), localWinSize=500, globalNoiseTh=0.75, 
 			smoothMethod=c('soft', 'hard'), method=c('dwt', 'modwt')) {

	if (!requireNamespace("waveslim", quietly = TRUE)) {
        stop('Please install the waveslim package to use smoothDWT()')
    }
	smoothMethod <- match.arg(smoothMethod)
	method <- match.arg(method)
	specLength <- length(ms)

	ms <- extendNBase(ms, nLevel=nLevel, method='open', direction='right')
	if (method == 'dwt') {
		coef <- waveslim::dwt(ms, wf="la8", n.levels=nLevel, boundary="reflection")		
	} else {
		coef <- waveslim::modwt(ms, wf="la8", n.levels=nLevel, boundary="reflection")		
	}
	
	localNoiseTh[localNoiseTh > 1] <- 1
	globalNoiseTh[globalNoiseTh > 1] <- 1
	len <- length(localNoiseTh)
	if (len < nLevel) {
		localNoiseTh <- c(localNoiseTh, rep(localNoiseTh[len], nLevel - len))
	}
	len <- length(globalNoiseTh)
	if (len < nLevel) {
		globalNoiseTh <- c(globalNoiseTh, rep(globalNoiseTh[len], nLevel - len))
	}
	
	coef.new <- coef
	for (i in 1:nLevel) {
		if (localNoiseTh[i] == 1 | globalNoiseTh[i] == 1 ) {
			coef.new[[i]][] <- 0
		} else {
			coef.i <- coef[[i]]
			## Two thresholds are used here. One is the global threshold, another one is local threshold
			## The local peaks should be above both thresholds
			globalTh.i <- stats::quantile(abs(coef.i), globalNoiseTh[i])
			ind.global <- which(abs(coef.i) > globalTh.i)

			## Transform the vector as a matrix with column length equals winSize
			##		and find the maximum position at each row.
			len <- length(coef.i)
			temp <- matrix(extendNBase(abs(coef.i), base=localWinSize, nLevel=1, method='open', direction='right'), nrow=localWinSize)

			localTh.i <- apply(temp, 2, function(x) stats::quantile(x, localNoiseTh[i]))
			localTh.i <- rep(1, localWinSize) %*% t(localTh.i)
			ind.local <- which(temp > localTh.i)
			selInd <- ind.local[ind.local %in% ind.global]

			coef.i[] <- 0
			if (smoothMethod == 'soft') {
				coef.i[selInd] <- sign(coef[[i]][selInd])*(abs(coef[[i]][selInd]) - localTh.i[selInd])
			} else {
				coef.i[selInd] <- coef[[i]][selInd]
			}
			coef.new[[i]] <- coef.i
		}
	}

	#coefMatrix <- matrix(unlist(coef.new), ncol=nLevel + 1)
	#approx <- coefMatrix[1:specLength, nLevel + 1]
	#detail <- rowSums(coefMatrix[1:specLength, 1:nLevel])
	temp <- coef.new
	temp[[nLevel + 1]][] <- 0
	if (method == 'dwt') {
		detail <- waveslim::idwt(temp)[1:specLength]		
	} else {
		detail <- waveslim::imodwt(temp)[1:specLength]		
	}
	temp <- coef.new
	for (i in 1:nLevel)  temp[[i]][] <- 0	
	if (method == 'dwt') {
		approx <- waveslim::idwt(temp)[1:specLength]		
		smoothMS <- waveslim::idwt(coef.new)[1:specLength]
	} else {
		approx <- waveslim::imodwt(temp)[1:specLength]		
		smoothMS <- waveslim::imodwt(coef.new)[1:specLength]
	}
	
	## Inverse modwt transform	
	#smoothMS <- smoothMS[1:specLength]
	attr(smoothMS, 'approximate') <- approx
	## Approximate can be calculated by: approx <- smoothMS - detail
	attr(smoothMS, 'detail') <- detail
	return(smoothMS)
}
