\name{getRidge}
\alias{getRidge}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Identify ridges based on the local maximum matrix}
\description{
Identify ridges by connecting the local maximum of 2-D CWT coefficients from the coarse scale to detail scale. The local maximum matrix is returned  from \code{\link{getLocalMaximumCWT}}
}
\usage{
getRidge(localMax, iInit = ncol(localMax), step = -1, iFinal = 1, minWinSize= 5, gapTh = 3, skip = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{localMax}{ The local maximum matrix is returned  from \code{\link{getLocalMaximumCWT}} with 1 represents maximum, others are 0. }
  \item{iInit}{ The start column to search ridge. By default, it starts from the coarsest scale level.}
  \item{step}{ Search step. -1 by default, which means searching from coarse scale to detail scale column by column. }
  \item{iFinal}{ The final column index of search ridge. }
  \item{minWinSize}{ The minimum slide window size used. }
  \item{gapTh}{ The gap allowed during searching for ridge. 3 by default. }
  \item{skip}{ The column to be skipped during search. }
}

\value{
 Return a list of ridge, with m/z index of peak at the raw spectrum as the names of the ridges.
}
\references{ 
Pan Du, Warren A. Kibbe, Simon M. Lin, Improved peak detection in mass spectrum by incorporating continuous wavelet transform-based pattern matching, Bioinformatics Advance Access published online on July 4, 2006
(http://bioinformatics.oxfordjournals.org/cgi/content/abstract/btl355?ijkey=HfCb7kHvQBY6HE6&keytype=ref)
}
\author{ Pan Du, Simon Lin }
\seealso{ \code{\link{getLocalMaximumCWT}}, \code{\link{identifyMajorPeaks}}}

\examples{
	data(exampleMS)
	scales <- seq(1, 64, 3)
	wCoefs <- cwt(exampleMS[5000:11000], scales=scales, wavelet='mexh')
	
	localMax <- getLocalMaximumCWT(wCoefs)
	ridgeList <- getRidge(localMax)
	plotRidgeList(ridgeList)
}
\keyword{methods}
