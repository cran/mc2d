\name{multinormal}
\alias{multinormal}
\alias{rmultinormal}
\title{
The Vectorized Multivariate Random Deviates
}
\description{
This function is the vectorized version of the \samp{rmvnorm} from
the \samp{mvtnorm} library. It provides a random number generator for
the multivariate normal distribution with varying vectors of means
and varying covariance matrixes.
}
\usage{
rmultinormal(n, mean, sigma, method=c("eigen", "svd", "chol"))
}
\examples{
(mean <- c(10, 0))
(sigma <- matrix(c(1, 2, 2, 10), ncol=2))
(sigma <- as.vector(sigma))
round(rmultinormal(10, mean, sigma)) 

(mean <- matrix(c(10, 0, 0, 10), ncol=2))
round(rmultinormal(10, mean, sigma))

(mean <- c(10, 0))
(sigma <- matrix(c(1, 2, 2, 10, 10, 2, 2, 1), nrow=2, byrow=TRUE))
round(rmultinormal(10, mean, sigma))

(mean <- matrix(c(10, 0, 0, 10), ncol=2))
(sigma <- matrix(c(1, 2, 2, 10, 10, 2, 2, 1), nrow=2, byrow=TRUE))
round(rmultinormal(10, mean, sigma))

(mean <- c(10, 0))
(sigma <- matrix(c(1, 2, 2, 10, 10, 2, 2, 1), nrow=2, byrow=TRUE))
round(x <- rmultinormal(1000, mean, sigma))
plot(x)
}
\arguments{
\item{n}{Number of observations.}
\item{mean}{Vector of means (if unique for all \samp{n}) or array of
means (if varying according to \samp{n}).}
\item{sigma}{Covariance vector corresponding to the coercion of the
covariance matrix into a vector (if unique for all \samp{n}) or array
of covariance vectors (if varying according to \samp{n}).}
\item{method}{Matrix decomposition used to determine the matrix root
of sigma, possible methods are eigenvalue decomposition ("eigen",
default), singular value decomposition ("svd"), and Cholesky
decomposition ("chol").}
}
\details{
\samp{rmvnorm(n, m, s)} is equivalent to \samp{rmultinormal(n, m,
as.vector(s))}.

If \samp{mean} and/or \samp{sigma} is a matrix, the first random
deviate will use the first row of \samp{mean} and/or \samp{sigma},
the second random deviate will use the second row of \samp{mean}
and/or \samp{sigma}, ... recycling being permitted by raw. If
\samp{mean} is a vector of length \samp{l} or is a matrix with
\samp{l} columns, \samp{sigma} should be a vector of length \samp{l x
l} or a matrix of number of \samp{l x 2} columns.
}
\note{
The use of a varying sigma may be very time consumming.
}
\keyword{distribution}
