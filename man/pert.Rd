\name{pert}
\alias{pert}
\alias{dpert}
\alias{ppert}
\alias{qpert}
\alias{rpert}
\title{
The Pert Distribution
}
\description{
Density, distribution function, quantile function and random
generation for the pert distribution with minimum equal to
\samp{min}, mode equal to \samp{mode} and maximum equal to
\samp{max}.
}
\usage{
dpert(x, min=-1, mode=0, max=1, shape=4, log=FALSE)
ppert(q, min=-1, mode=0, max=1, shape=4, lower.tail=TRUE, log.p=FALSE)
qpert(p, min=-1, mode=0, max=1, shape=4, lower.tail=TRUE, log.p=FALSE)
rpert(n, min=-1, mode=0, max=1, shape=4)
}
\examples{
curve(dpert(x, min=3, mode=5, max=10, shape=6), from = 2, to = 11, lty=3)
curve(dpert(x, min=3, mode=5, max=10), from = 2, to = 11, add=TRUE)
curve(dpert(x, min=3, mode=5, max=10, shape=2), from = 2, to = 11, add=TRUE, lty=2)
legend(x = 8, y = 2, c("Default", "shape:2", "shape:6"), lty=1:3)



}
\arguments{
\item{x,q}{Vector of quantiles.}
\item{p}{Vector of probabilities.}
\item{n}{Number of observations. If length(n) \eqn{>} 1, the length
is taken to be the number required.}
\item{min}{Vector of minima.}
\item{mode}{Vector of modes.}
\item{max}{Vector of maxima.}
\item{shape}{Vector of scaling parameters.}
\item{log, log.p}{Logical; if \samp{TRUE}, probabilities \samp{p} are
given as \samp{log(p)}.}
\item{lower.tail}{Logical; if \samp{TRUE} (default), probabilities
are \samp{P[X <= x]}, otherwise, \samp{P[X > x]}.}
}
\details{
The Pert distribution is a special case of the Beta distribution
specified by the following parameters. Given:
\deqn{\mu=\frac{min+max+shape\times mode}{shape+2}}{mu = (min + max +
shape * mode)/(shape + 2)} the values of \eqn{\alpha_{1}}{shape1} and
\eqn{\alpha_{2}}{shape2} are \deqn{\alpha_{1}=\frac{(\mu-min)(2\times
mode-min-max)}{(mode-\mu)(max-min)}}{shape1=(mu - min)*(2
mode-min-max)/((mode-mu)*(max - min)}

\deqn{\alpha_{2}=\frac{\alpha_{1}\times
(max-\mu)}{mu-min}}{shape2=shape1*(max - mu)/(mu - min)}

on the domain \samp{[min, max]}.

 If \eqn{\mu=mode}{mu=mode}, \eqn{\alpha_{1}}{shape1} is set to
\eqn{1+\nu/2}{1+shape/2}.
}
\value{
\samp{dpert} gives the density, \samp{ppert} gives the distribution
function, \samp{qpert} gives the quantile function, and \samp{rpert}
generates random deviates.
}
\references{
Vose D. Risk Analysis - A Quantitative Guide (John Wiley \& Sons,
2000).
}
\author{
Regis Pouillot
}
\seealso{
\code{\link{Beta}}
}
\keyword{distribution}
