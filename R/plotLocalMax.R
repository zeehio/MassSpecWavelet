"plotLocalMax" <-
function(localMax, wCoefs=NULL, range=c(1, nrow(localMax)), colorMap='RYB', main=NULL, cex=3, pch='.', ... ) {
	if (colorMap == 'RYB') {
		rgb.palette <- grDevices::colorRampPalette(c("red", "orange", "blue"), space = "rgb")
		colorMap <- rgb.palette(255)
		#colorMap <- heat.colors(255)
	}
	
	## Check range parameter
	range <- sort(range, decreasing=FALSE)
	if (length(range) == 1) range <- c(range, ncol(localMax))
	if (range[1] < 1) range[1] <- 1
	if (range[2] > nrow(localMax)) range[2] <- ncol(localMax)

	test <- localMax[range[1]:range[2],]
	tR <- range[1] - 1 + row(test)
	tC <- col(test)
	ind <- which(test > 0)
	if (is.null(wCoefs)) {
		col <- 'blue'
	} else {
		selRidgeValue <- wCoefs[range[1]:range[2], ][ind]
		maxSelV <- log(1 + max(selRidgeValue))
		col <- colorMap[1+log(1+abs(selRidgeValue))*255/maxSelV]		
	}

	plot(seq(range[1],range[2], length=ncol(localMax)), 1:ncol(localMax), type='n', xlab='m/z index', ylab='CWT coefficient scale', main=main)
	ind <- which(localMax[range[1]:range[2],] > 0)
	graphics::points(tR[ind], tC[ind], col=col, cex=cex, pch=pch, ...)
}

