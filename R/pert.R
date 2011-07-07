#<<BEGIN>>
dpert <- function(x,min=-1,mode=0,max=1,shape=4,log=FALSE)
#TITLE The Pert Distribution
#NAME pert
#KEYWORDS distribution
#DESCRIPTION
#Density, distribution function, quantile function and random generation
#for the pert distribution with minimum equal to \samp{min}, mode equal to \samp{mode}
#and maximum equal to \samp{max}.
#INPUTS
#{x,q}<<Vector of quantiles.>>
#{p}<<Vector of probabilities.>>
#{n}<<Number of observations. If length(n) > 1, the length is taken to be the number required.>>
#[INPUTS]
#{min}<<Vector of minima.>>
#{mode}<<Vector of modes.>>
#{max}<<Vector of maxima.>>
#{shape}<<Vector of scaling parameters.>>
#{log, log.p}<<Logical; if \samp{TRUE}, probabilities \samp{p} are given as \samp{log(p)}.>>
#{lower.tail}<<Logical; if \samp{TRUE} (default), probabilities are \samp{P[X <= x]}, otherwise, \samp{P[X > x]}.>>
#DETAILS
#The Pert distribution is a special case of the Beta distribution specified by the following parameters.
#Given:
#\deqn{\mu=\frac{min+max+shape\times mode}{shape+2}}{mu = (min + max + shape * mode)/(shape + 2)}
#the values of \eqn{\alpha_{1}}{shape1} and \eqn{\alpha_{2}}{shape2} are
#\deqn{\alpha_{1}=\frac{(\mu-min)(2\times mode-min-max)}{(mode-\mu)(max-min)}}{shape1=(mu - min)*(2 mode-min-max)/((mode-mu)*(max - min)}
#
#\deqn{\alpha_{2}=\frac{\alpha_{1}\times (max-\mu)}{mu-min}}{shape2=shape1*(max - mu)/(mu - min)}
#
#on the domain \samp{[min, max]}.</>
#If \eqn{\mu=mode}{mu=mode}, \eqn{\alpha_{1}}{shape1} is set to \eqn{1+\nu/2}{1+shape/2}.
#REFERENCE
#Vose D. Risk Analysis - A Quantitative Guide (John Wiley & Sons, 2000).
#EXAMPLE
#curve(dpert(x,min=3,mode=5,max=10,shape=6), from = 2, to = 11, lty=3)
#curve(dpert(x,min=3,mode=5,max=10), from = 2, to = 11, add=TRUE)
#curve(dpert(x,min=3,mode=5,max=10,shape=2), from = 2, to = 11, add=TRUE,lty=2)
#legend(x = 8, y = 2, c("Default","shape:2","shape:6"), lty=1:3)
#SEE ALSO
#\code{\link{Beta}}
#VALUE
#\samp{dpert} gives the density, \samp{ppert} gives the distribution function,
#\samp{qpert} gives the quantile function, and \samp{rpert} generates random deviates.

#CREATED 08-02-20
#--------------------------------------------
{
	if(length(x) == 0) return(numeric(0))
	
	mu <- (min+max+shape*mode)/(shape+2)
	a1 <- ifelse(mapply(function(x,y) isTRUE(all.equal(x,y)),mu,mode),
                1+shape/2,
                (mu-min)*(2*mode-min-max)/((mode-mu)*(max-min)))
	a2 <- a1*(max-mu)/(mu-min)
	oldw <- options(warn = -1)
	
	d <- dbeta(x=(x-min)/(max-min),shape1=a1,shape2=a2) / (max-min)
	options(warn = oldw$warn)
	
	d[x < min | x > max] <- 0
	d[mode < min | max < mode] <- NaN
	d[shape <= -2] <- NaN
	if(log) d <- log(d)
	if(any(is.na(d))) warning("NaN in dpert")
  return(d)}

#<<BEGIN>>
ppert <- function(q,min=-1,mode=0,max=1,shape=4,lower.tail = TRUE, log.p = FALSE)
#ISALIAS dpert
#--------------------------------------------
{
	if(length(q) == 0) return(numeric(0))
	
	mu <- (min + max + shape*mode)/(shape + 2)
	a1 <- ifelse(mapply(function(x,y) isTRUE(all.equal(x,y)),mu,mode),
                1+shape/2,
                (mu-min)*(2*mode-min-max)/((mode-mu)*(max-min)))
	a2 <- a1*(max-mu)/(mu-min)
	oldw <- options(warn = -1)
	p <- pbeta(q=(q-min)/(max-min),shape1=a1,shape2=a2)
	options(warn = oldw$warn)
	p[q < min] <- 0
	p[q >= max] <- 1
	p[mode < min | max < mode] <- NaN
	p[shape <= -2] <- NaN
	if(!lower.tail) p <- 1-p
	if(log.p) p <- log(p)
	if(any(is.na(p))) warning("NaN in ppert")
  return(p)}

#<<BEGIN>>
qpert <- function(p,min=-1,mode=0,max=1,shape=4,lower.tail=TRUE,log.p=FALSE)
#ISALIAS dpert
#--------------------------------------------
{
  if(length(p) == 0) return(numeric(0))
  if(log.p) p <- exp(p)
  if(!lower.tail) p <- 1 - p
	mu <- (min+max+shape*mode)/(shape+2)
	a1 <- ifelse(mapply(function(x,y) isTRUE(all.equal(x,y)),mu,mode),
                1+shape/2,
                (mu-min)*(2*mode-min-max)/((mode-mu)*(max-min)))
	a2 <- a1*(max-mu)/(mu-min)
	oldw <- options(warn = -1)
  q <- qbeta(p, shape1=a1, shape2=a2)
  options(warn = oldw$warn)
  q <- q * (max-min) + min
  #Very special case min = max = mode
  q[mapply(function(x,y) isTRUE(all.equal(x,y)),min,max)] <- 1
  q[p < 0 | p > 1] <- NaN
  q[mode < min | max < mode] <- NaN
  q[shape <= -2] <- NaN
  if(any(is.na(q))) warning("NaN in qpert")
  return(q)}


#<<BEGIN>>
rpert <- function(n,min=-1,mode=0,max=1,shape=4)
#ISALIAS dpert
#--------------------------------------------
{ if(length(n) > 1) n <- length(n)
  if(length(n) == 0 || as.integer(n) == 0) return(numeric(0))
  n <- as.integer(n)
  if(n < 0) stop("integer(n) can not be negative in rpert")
  
  oldw <- options(warn = -1)
  r <- qpert(runif(n),min=min,mode=mode,max=max,shape=shape,lower.tail=TRUE,log.p=FALSE)
  options(warn = oldw$warn)
  if(any(is.na(r))) warning("NaN in rpert")
  return(r)
}
