\name{triangular}
\alias{triangular}
\alias{dtriang}
\alias{ptriang}
\alias{qtriang}
\alias{rtriang}
\title{
The Triangular Distribution
}
\description{
Density, distribution function, quantile function and random
generation for the triangular distribution with minimum equal to
\samp{min}, mode equal \samp{mode} and maximum equal to \samp{max}.
}
\usage{
dtriang(x, min=-1, mode=0, max=1, log=FALSE)
ptriang(q, min=-1, mode=0, max=1, lower.tail=TRUE, log.p=FALSE)
qtriang(p, min=-1, mode=0, max=1, lower.tail=TRUE, log.p=FALSE)
rtriang(n, min=-1, mode=0, max=1)
}
\examples{
curve(dtriang(x, min=3, mode=5, max=10), from = 2, to = 11)



}
\arguments{
\item{x,q}{vector of quantiles.}
\item{p}{vector of probabilities.}
\item{n}{number of observations. If length(n) \eqn{>} 1, the length
is taken to be the number required.}
\item{min}{vector of minima.}
\item{mode}{vector of modes.}
\item{max}{vector of maxima.}
\item{log, log.p}{logical; if \samp{TRUE}, probabilities \samp{p} are
given as \samp{log(p)}.}
\item{lower.tail}{logical; if \samp{TRUE} (default), probabilities
are \samp{P[X <= x]}, otherwise, \samp{P[X > x]}.}
}
\value{
\samp{dtriang} gives the density, \samp{ptriang} gives the
distribution function, \samp{qtriang} gives the quantile function,
and \samp{rtriang} generates random deviates.
}
\author{
Regis Pouillot
}
\keyword{distribution}
