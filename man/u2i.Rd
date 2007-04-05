\name{u2i}
\alias{u2i}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Transfer m/z value to m/z index}
\description{
  Transfer m/z value to m/z index
}
\usage{
u2i(u, t0 = 9e-08, a = 336301655.051424, b = 0, U = 20000, delta = 4e-09)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{u}{ a vector of m/z values }
  \item{t0}{ initial time }
  \item{a}{ accelerate rate }
  \item{b}{ constant shift }
  \item{U}{ Unit }
  \item{delta}{sample interval in second }
}
\details{
This function is internally used.
It is used when the data has not been calibrated or no m/z value information available.
The default parameters are just for the example data set. Different machines have different parameters.

Here is the conversion equation:
 \deqn{i = \frac{(\frac{u/U-b}{a})^{1/2} + t0}{delta}}{i = (((u/U - b)/a)^(1/2) + t0)/delta}
}
\value{
 a vector of m/z indexes
}
\author{ Simon Lin }

\seealso{ \code{\link{i2u}} }

\keyword{methods}
