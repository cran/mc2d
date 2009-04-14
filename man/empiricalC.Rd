\name{empiricalC}
\alias{empiricalC}
\alias{dempiricalC}
\alias{pempiricalC}
\alias{qempiricalC}
\alias{rempiricalC}
\title{
The Continuous Empirical Distribution
}
\description{
Density, distribution function and random generation for a continuous
empirical distribution.
}
\usage{
dempiricalC(x, min, max, values, prob=rep(1, length(values)), log=FALSE)
pempiricalC(q, min, max, values, prob=rep(1, length(values)), lower.tail=TRUE,
	  log.p=FALSE)
qempiricalC(p, min, max, values, prob=rep(1, length(values)), lower.tail=TRUE,
	  log.p=FALSE)
rempiricalC(n, min, max, values, prob=rep(1, length(values)))
}
\examples{
prob <- c(2, 3, 1, 6, 1)
values <- 1:5
par(mfrow=c(1, 2))
curve(dempiricalC(x, min=0, max=6, values, prob), from=-1, to=7, n=1001)
curve(pempiricalC(x, min=0, max=6, values, prob), from=-1, to=7, n=1001)



}
\arguments{
\item{x, q}{Vector of quantiles.}
\item{p}{Vector of probabilities.}
\item{n}{Number of random values. If \samp{length(n) > 1}, the length
is taken to be the number required.}
\item{min}{A finite minimal value.}
\item{max}{A finite maximal value.}
\item{values}{Vector of numerical values.}
\item{prob}{Optionnal vector of count or probabilities.}
\item{log, log.p}{logical; if \samp{TRUE}, probabilities \samp{p} are
given as \samp{log(p)}.}
\item{lower.tail}{logical; if \samp{TRUE} (default), probabilities
are \samp{P[X <= x]}, otherwise, \samp{P[X > x]}.}
}
\details{
Given \eqn{p_{i}}{p_i}, the distribution value for \eqn{x_{i}}{x_i}
with \samp{i} the rank \eqn{i = 0, 1, 2, \ldots, N+1},
\eqn{x_{0}=min}{x_0 = min} and \eqn{x_{N+1}=max}{x_(N+1) = max} the
density is:
\deqn{f(x)=p_{i}+(\frac{x-x_{i}}{x_{i+1}-x_{i}})(p_{i+1}-p_{i})}{f(x)
= p_i + (p_(i+1) - p_i)/(x_(i+1) - x_i) for x_i<=x<x_(i+1)} The
\samp{p} values being normalized to give the distribution a unit
area.

This function is not vectorized. It can not currently be use with
varying parameters.
}
\value{
\samp{dempiricalC} gives the density, \samp{pempiricalC} gives the
distribution function, \samp{qempiricalC} gives the quantile function
and \samp{rempiricalC} generates random deviates.
}
\author{
Regis Pouillot
}
\seealso{
\code{\link{empiricalD}}
}
\keyword{distribution}
