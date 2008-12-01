#<<BEGIN>>
dmultinomial <- function (x, size = NULL, prob, log = FALSE)
#TITLE The Vectorized Multinomial Distribution
#KEYWORDS distribution
#DESCRIPTION
#Generate multinomially distributed random number vectors and compute multinomial probabilities.
#INPUTS
#{x}<<Vector of length K of integers in 0:size.>>
#{n}<<Number of random vectors to draw.>>
#{size}<<A vector of integers, say N, specifying the total number of objects that are put
#into K boxes in the typical multinomial experiment. For \samp{dmultinom}, it defaults to \samp{sum(x)}.
#The first element correspond to the vector \samp{prob} or the first row of \samp{prob}, ...>>
#{prob}<<Numeric non-negative vector of length K, or matrix of size \samp{(x x K)}
#specifying the probability for the K classes; is internally normalized to sum 1.>>
#{log}<<Logical; if TRUE, log probabilities are computed.>>
#EXAMPLE
#prob <- c(1,2,7)
#rmultinomial(4,1000,prob)
#rmultinomial(4,c(10,100,1000,10000),prob)
#
### rmultinomial used with mcstoc
### (uncertain size and prob)
#s <- mcstoc(rpois,"U",lambda=50)
#p <- mcstoc(rdirichlet,"U",nvariates=3,alpha=c(4,10,20))
#mcstoc(rmultinomial,"VU",nvariates=3,size=s, prob=p)
#DETAILS
#This function is the vectorized version of \code{\link{rmultinom}} and \code{\link{dmultinom}}.
#Recycling is permitted.
#--------------------------------------------
{
    if(is.vector(prob)) prob <- t(prob)
    if(is.vector(x)) prob <- t(x)
    K <- ncol(prob)
    n <- nrow(x)
    if (ncol(x) != K)
        stop("x[] and prob[] must be equal length vectors or equal col matrix.")
    if(nrow(prob)!= n) prob <- matrix(t(prob),ncol=K,nrow=n, byrow=TRUE)        #recycling
    if(!is.null(size) && nrow(size)!= n) size <- matrix(t(size),ncol=K,nrow=n, byrow=TRUE)        #recycling
    if (any(prob < 0) || any(s <- apply(prob,1,sum) == 0))
        stop("probabilities cannot be negative nor all 0.")
    prob <- prob/s
    x <- as.integer(x + 0.5)
    if (any(x < 0))
        stop("'x' must be non-negative")
    N <- apply(x,1,sum)
    if (is.null(size))
        size <- N
    else if (size != N)
        stop("size != sum(x), i.e. one is wrong")
    i0 <- prob == 0
    if (any(i0)) {
        if (any(x[i0] != 0))
            return(if (log) -Inf else 0)
        if (all(i0))
            return(if (log) 0 else 1)
        x <- x[!i0]
        prob <- prob[!i0]
    }

    r <- sapply(1:n,function(y) lgamma(size[y,] + 1) + sum(x[y,] * log(prob[y,]) - lgamma(x[y,] + 1)))
    if (log)
        r
    else exp(r)
}

#<<BEGIN>>
rmultinomial <- function(n, size, prob)
#ISALIAS dmultinomial
#--------------------------------------------
{
  if(length(n)!=1) n <- length(n)
  if(is.vector(prob) || (dim(prob)[1]) == 1) {
    if(length(size)==1) return(t(rmultinom(n,size,prob)))      # classical
    prob <- matrix(prob,nrow=1)
    }

  nrp <- nrow(prob)
  mnr <- min( max(nrp ,length(size)), n)
  ss  <- rep(size,length.out=mnr)              # recycling size
  
  if(nrp != mnr) prob <- matrix(t(prob),ncol=ncol(prob),nrow=mnr,byrow=TRUE)    # recycling prob

  n1 <- n%/%mnr
  n2 <- n%%mnr

  res <- sapply(1:mnr,function(x) rmultinom(n1,ss[x],prob[x,]))
  res <- matrix(res,ncol=ncol(prob),byrow=TRUE)
  index <- as.vector(matrix(1:(mnr*n1),ncol=mnr,byrow=TRUE))
  res <- res[index,]

  if (n2 != 0){
    res <- rbind(res,t(sapply(1:n2,function(x) rmultinom(1,ss[x],prob[x,]))))
    }
  return(res)
}

