\name{i2u}
\alias{i2u}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Transfer m/z index to m/z value }
\description{
  Transfer m/z index to m/z value
}
\usage{
i2u(i, t0 = 9e-08, a = 336301655.051424, b = 0, U = 20000, delta = 4e-09)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{i}{  a vector of m/z index }
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
 	t = i * delta, (unit in seconds),

 	\deqn{u = (a*(t-t0)^{2} + b)*U}{u = (a*(t-t0)^2 + b)*U}, (here unit is m/z at 150KDa)
}
\value{
 a vector of m/z values
}

\author{ Simon Lin }

\seealso{ \code{\link{u2i}} }

\keyword{methods}
