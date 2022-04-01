\name{sav.gol}
\alias{sav.gol}

\title{ Estimate the baseline by using Savitzky-Golay Algorithm }
\description{
  Estimate the baseline by using Savitzky-Golay Algorithm
}
\usage{
sav.gol(T, fl, forder = 4, dorder = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{T}{ vector of signals to be filtered }
  \item{fl}{ filter length (for instance fl = 51..151) }
  \item{forder}{ filter order (2 = quadratic filter, 4= quartic) }
  \item{dorder}{ derivative order (0 = smoothing, 1 = first derivative, etc.) }
}

\value{
The return is a smoothed vector (baseline).
}

\author{Steffen Neumann <sneumann@ipb-halle.de> }
\note{ 
This function was added by Steffen Neumann. We appreciated his help to make the package better.
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
