\name{cwt}
\alias{cwt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Continuous Wavelet Transform (CWT) }
\description{
  CWT(Continuous Wavelet Transform) with Mexican Hat wavelet (by default) to match the peaks in Mass Spectrometry spectrum
}
\usage{
cwt(ms, scales = 1, wavelet = "mexh")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ms}{ Mass Spectrometry spectrum (a vector of MS intensities)}
  \item{scales}{ a vector represents the scales at which to perform CWT. See the Details section}
  \item{wavelet}{The wavelet base, Mexican Hat by default. User can provide wavelet Psi(x) as a form of two row matrix. The first 
	row is the x value, and the second row is Psi(x) corresponding to x.  }
}

\value{
The return is the 2-D CWT coefficient matrix, with column names as the scale. Each column is the CWT coefficients at that scale.
}

\details{
The default mother wavelet is a Mexican Hat evaluated in the [-8,8] range using 1024 points (a step of roughly 1/64):

\deqn{ \psi(x) = \frac{2}{\sqrt{3}} \pi^{-0.25} ( 1 - x^2 ) \exp{-x^2/2} }

The \eqn{\sigma} of the mother Mexican Hat is of one x unit.

The scaled wavelet is a downsampled version of the mother wavelet. The \code{scale} determines how many samples per \eqn{x} unit
are taken. For instance, with the default Mexican Hat wavelet, a \code{scales = 2} will evaluate the [-8, 8] range sampling
twice per \eqn{x} unit, this is with a step of 0.5. This generates a 33 points long scaled wavelet. Choosing this type of scaling is convenient because the
scaled wavelet becomes a wavelet of \eqn{\sigma = \code{scales}} points. Using the default wavelet, a \code{scales} smaller than
 1 will show sampling issues, while a \code{scales} larger than 64 will resample points from the original mother wavelet. If you
 need to use \code{scales} larger than 64, consider providing your own mother wavelet. See the examples.


According to \doi{10.1063/1.3505103}, if your spectrum has a gaussian peak shape of variance \eqn{s^2} points then the optimal scale to detect it
should be \code{scales = 1.9 s}. If your spectrum has a Lorentzian peak shape of half-width-half-maximum \eqn{L} points then the
optimal scale should be \code{scales = 2.9 L}. Take these values as suggestions for your analysis.

}

\author{Pan Du, Simon Lin}

\examples{
	data(exampleMS)
	scales <- seq(1, 64, 3)
	wCoefs <- cwt(exampleMS[5000:11000], scales=scales, wavelet='mexh')

	## Plot the 2-D CWT coefficients as image (It may take a while!)
	xTickInterval <- 1000
	image(5000:11000, scales, wCoefs, col=terrain.colors(256), axes=FALSE,
          xlab='m/z index', ylab='CWT coefficient scale', main='CWT coefficients')
	axis(1, at=seq(5000, 11000, by=xTickInterval))
	axis(2, at=c(1, seq(10, 64, by=10)))
	box()

	## Provide a larger wavelet:
	scales <- c(seq(1, 64, 3), seq(72, 128, 8))
	psi_xval <- seq(-8, 8, length.out = 4096)
	psi <- (2/sqrt(3) * pi^(-0.25)) * (1 - psi_xval^2) *exp(-psi_xval^2/2)
	wCoefs <- cwt(exampleMS[5000:11000], scales=scales, wavelet=rbind(psi_xval, psi))
	xTickInterval <- 1000
	image(5000:11000, scales, wCoefs, col=terrain.colors(256), axes=FALSE,
          xlab='m/z index', ylab='CWT coefficient scale', main='CWT coefficients')
	axis(1, at=seq(5000, 11000, by=xTickInterval))
	axis(2, at=c(1, seq(10, 128, by=10)))
	box()

}
\keyword{ methods }
