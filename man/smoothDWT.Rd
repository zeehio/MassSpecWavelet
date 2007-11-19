\name{smoothDWT}
\alias{smoothDWT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ smooth (denoise) the spectrum by DWT (Discrete Wavelet Transform) }
\description{
  Smooth (denoise) the spectrum by DWT (Discrete Wavelet Transform)
}
\usage{
smoothDWT(ms, nLevel = 6, wf = "la8", localNoiseTh = seq(1, 0, by = -0.2), localWinSize = 500, globalNoiseTh = 0.75, smoothMethod = c("soft", "hard"), method = c('dwt', 'modwt'))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ms}{ a vector representing the mass spectrum }
  \item{nLevel}{ the level of DWT decomposition }
  \item{wf}{ the name of wavelet for DWT }
  \item{localNoiseTh}{ local noise level threshold }
  \item{localWinSize}{ local window size for estimate local noise threshold }
  \item{globalNoiseTh}{ global noise level threshold }
  \item{smoothMethod}{ the method used for denoising. 'hard' means keeping the dwt coefficients higher than the threshold unchanged; "soft" means the dwt coefficients higher than the threshold were subtracted by the threshold. }
  \item{method}{ 'dwt' or 'modwt' used for decomposition }
}
\details{
  

}
\value{
return the smoothed mass spectrum with the 'detail' component of DWT as an attribute 'detail'.
}
\references{  }
\author{ Pan Du }

\examples{

}
\keyword{ methods }

