% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/betasubjective.R
\encoding{UTF-8}
\name{BetaSubjective}
\alias{BetaSubjective}
\alias{dbetasubj}
\alias{pbetasubj}
\alias{qbetasubj}
\alias{rbetasubj}
\title{The BetaSubjective Distribution}
\usage{
dbetasubj(x, 
  min,
  mode,
  mean,
  max, 
  log = FALSE)

pbetasubj(q, 
  min,
  mode,
  mean,
  max, 
  lower.tail = TRUE,
  log.p = FALSE
)

qbetasubj(p, 
  min,
  mode,
  mean,
  max, 
  lower.tail = TRUE, 
  log.p = FALSE
)

rbetasubj(n, 
  min,
  mode,
  mean,
  max
)

pbetasubj(q, min, mode, mean, max, lower.tail = TRUE, log.p = FALSE)

qbetasubj(p, min, mode, mean, max, lower.tail = TRUE, log.p = FALSE)

rbetasubj(n, min, mode, mean, max)
}
\arguments{
\item{x, q}{Vector of quantiles.}

\item{min}{continuous boundary parameter min < max}

\item{mode}{continuous parameter \eqn{min < mode < max} and \eqn{mode \ne mean}.}

\item{mean}{continuous parameter min < mean < max}

\item{max}{continuous boundary parameter}

\item{log, log.p}{Logical; if TRUE, probabilities p are given as log(p).}

\item{lower.tail}{Logical; if TRUE (default), probabilities are \eqn{P[X \le x]} otherwise, \eqn{P[X > x]}.}

\item{p}{Vector of probabilities.}

\item{n}{Number of observations.}
}
\description{
Density, distribution function, quantile function and random generation
for the "Beta Subjective" distribution
}
\details{
The Subjective beta distribution specifies a [stats::dbeta()] distribution defined by the minimum, most likely (mode), mean
and maximum values and can be used for fitting data for a variable that is bounded to the interval \eqn{[min, max]}. 
The shape parameters are calculated from the mode value and mean parameters. It can also be used to represent 
uncertainty in subjective expert estimates.

Define
\deqn{mid=(min+max)/2}{mid=(min+max)/2}
\deqn{a_{1}=2*\frac{(mean-min)*(mid-mode)}{((mean-mode)*(max-min))}}{a_1=2*(mean-min)*(mid-mode)/((mean-mode)*(max-min))}
\deqn{a_{2}=a_{1}*\frac{(max-mean)}{(mean-min)}}{a_2=a_1*(max-mean)/(mean-min)}

The subject beta distribution is a [stats::dbeta()] distribution defined on the \eqn{[min, max]} domain 
with parameter \eqn{shape1 = a_{1}} and \eqn{shape2 = a_{2}}.  

# Hence, it has density
# \deqn{f(x)=(x-min)^{(a_{1}-1)}*(max-x)^{(a_{2}-1)} / (B(a_{1},a_{2})*(max-min)^{(a_{1}+a_{2}-1)})}{f(x)=(x-min)^(a_1-1)*(max-x)^(a_2-1)/(B(a_1,a_2)*(max-min)^(a_1+a_2-1))}

# The cumulative distribution function is
# \deqn{F(x)=B_{z}(a_{1},a_{2})/B(a_{1},a_{2})=I_{z}(a_{1},a_{2})}{F(x)=B_z(a_1,a_2)/B(a_1,a_2)=I_z(a_1,a_2)}
# where \eqn{z=(x-min)/(max-min)}. Here B is the beta function and \eqn{B_z} is the incomplete beta function.

The parameter restrictions are:
\deqn{min <= mode <= max}{min <= mode <= max}
\deqn{min <= mean <= max}{min <= mean <= max}
If \eqn{mode > mean} then \eqn{mode > mid}, else \eqn{mode < mid}.
}
\examples{
curve(dbetasubj(x, min=0, mode=1, mean=2, max=5), from=-1,to=6) 
pbetasubj(q = seq(0,5,0.01), 0, 1, 2, 5)
qbetasubj(p = seq(0,1,0.01), 0, 1, 2, 5)
rbetasubj(n = 1e7, 0, 1, 2, 5)
}
\author{
Yu Chen
}
\keyword{distribution}
