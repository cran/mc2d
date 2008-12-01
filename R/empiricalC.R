#<<BEGIN>>
dempiricalC <- function(x, min, max, values, prob=rep(1,length(values)), log=FALSE)
#TITLE The Continuous Empirical Distribution
#NAME empiricalC
#KEYWORDS distribution
#DESCRIPTION
#Density, distribution function and random generation for a continuous empirical distribution.
#INPUTS
#{x, q}<<Vector of quantiles.>>
#{p}<<Vector of probabilities.>>
#{n}<<Number of random values. If \samp{length(n) > 1}, the length is taken to be the number required.>>
#{min}<<A finite minimal value.>>
#{max}<<A finite maximal value.>>
#{values}<<Vector of numerical values.>>
#[INPUTS]
#{prob}<<Optionnal vector of count or probabilities.>>
#{log, log.p}<<logical; if \samp{TRUE}, probabilities \samp{p} are given as \samp{log(p)}.>>
#{lower.tail}<<logical; if \samp{TRUE} (default), probabilities are \samp{P[X <= x]}, otherwise, \samp{P[X > x]}.>>
#DETAILS
#Given \eqn{p_{i}}{p_i}, the distribution value for \eqn{x_{i}}{x_i} with \samp{i} the rank \eqn{i = 0, 1, 2, \ldots, N+1},
# \eqn{x_{0}=min}{x_0 = min} and \eqn{x_{N+1}=max}{x_(N+1) = max} the density is:
# \deqn{f(x)=p_{i}+(\frac{x-x_{i}}{x_{i+1}-x_{i}})(p_{i+1}-p_{i})}{f(x) = p_i + (p_(i+1) - p_i)/(x_(i+1) - x_i) for x_i<=x<x_(i+1)}
# The \samp{p} values being normalized to give the distribution a unit area.
#
#This function is not vectorized. It can not currently be use with varying parameters.
#SEE ALSO
#\code{\link{empiricalD}}
#VALUE
#\samp{dempiricalC} gives the density, \samp{pempiricalC} gives the distribution function,
#\samp{qempiricalC} gives the quantile function and \samp{rempiricalC} generates random deviates.
#EXAMPLE
#prob <- c(2,3,1,6,1)
#values <- 1:5
#par(mfrow=c(1,2))
#curve(dempiricalC(x, min=0, max=6, values, prob), from=-1, to=7, n=1001)
#curve(pempiricalC(x, min=0, max=6, values, prob), from=-1, to=7, n=1001)
#AUTHOR Regis Pouillot
#CREATED 08-02-20
#--------------------------------------------
{
  
  
  
  if(min > max | min > min(values) | max < max(values) | !is.finite(min) | !is.finite(max)) stop("Error in min or max")

  val2 <- sort(unique(values))
  probi <- dempiricalD(val2, values=values, prob=prob, log=FALSE)

  prob <- c(0,0,probi,0,0)
  val <- c(-Inf,min,val2,max,Inf)
  
  h <- c(val2,max)-c(min,val2)
  a <- c(0,probi)
  b <- c(probi,0)
  Integ <- sum(h*(a+b)/2)

  d <- rep(NA,length(x))
  d[x > max | x < min] <- 0
  lesquel <- which(!is.na(x) & x <= max & x >= min)
  quel <- sapply(x[lesquel],function(y) min(which(y < val)))
  d[lesquel] <- prob[quel-1]+(x[lesquel]-val[quel-1])/(val[quel]-val[quel-1])*(prob[quel]-prob[quel-1])

  if(log) d <- log(d)
	if(any(is.na(d))) warning("NaN in dempiricalC")
  return(d/Integ)}

#<<BEGIN>>
pempiricalC <- function(q, min, max, values, prob=rep(1,length(values)),lower.tail = TRUE, log.p = FALSE)                 
#ISALIAS dempiricalC
#--------------------------------------------
{
  if(min > max | min > min(values) | max < max(values) | !is.finite(min) | !is.finite(max)) stop("Error in min or max")

  val2 <- sort(unique(values))
  probi   <- dempiricalC(val2,min=min,max=max,values=values,prob=prob,log=FALSE)

  h <- c(val2,max)-c(min,val2)
  a <- c(0,probi)
  b <- c(probi,0)
  probcum <- cumsum(h*(a+b)/2)

  probi <- c(0,0,probi,0,0)
  probcum <- c(0,0,probcum,1)
  val <- c(-Inf,min,val2,max,Inf)

  p <- rep(NA,length(q))
  p[q >= max] <- 1
  p[q <= min] <- 0
  lesquel <- which(!is.na(q) & q < max & q > min)
  quel <- sapply(q[lesquel],function(y) min(which(y < val)))

  p[lesquel] <- probcum[quel-1]+(q[lesquel]-val[quel-1])*
                (probi[quel-1]+((probi[quel]-probi[quel-1])*(q[lesquel]-val[quel-1])/(2*(val[quel]-val[quel-1]))))

  if(!lower.tail) p <- 1-p
  if(log.p) p <- log(p)
	if(any(is.na(p))) warning("NaN in pempiricalC")
  return(p)}


#<<BEGIN>>
qempiricalC <- function(p, min, max, values, prob=rep(1,length(values)), lower.tail = TRUE, log.p = FALSE)
#ISALIAS dempiricalC
#--------------------------------------------
{
  if(log.p) p <- exp(p)
  if(!lower.tail) p <- 1-p
  if(min > max | min > min(values) | max < max(values) | !is.finite(min) | !is.finite(max)) stop("Error in min or max")

  val2 <- sort(unique(values))
  probi   <- dempiricalC(val2,min=min,max=max,values=values,prob=prob,log=FALSE)

  h <- c(val2,max)-c(min,val2)
  a <- c(0,probi)
  b <- c(probi,0)
  probcum <- cumsum(h*(a+b)/2)

  probi <- c(0,0,probi,0,0)
  probcum <- c(0,0,probcum,Inf)
  val <- c(-Inf,min,val2,max,Inf)

  q <- rep(NA,length(p))
  q[p > 1] <- NaN
  q[p == 1] <- max
  q[p < 0] <- NaN

  lesquel <- which(!is.na(p) & p < 1 & p >= 0)
  
  quel <- sapply(p[lesquel],function(y) min(which(y < probcum)))

  a <- (probi[quel]-probi[quel-1])/(val[quel]-val[quel-1])/2
  b <- probi[quel-1]
  c <- probcum[quel-1]-p[lesquel]
  d <- b^2-4*a*c
  q[lesquel] <- (-b+sqrt(d))/2/a  + val[quel-1]
  
	if(any(is.na(q))) warning("NaN in qempiricalC")
  return(q)}

#<<BEGIN>>
rempiricalC <- function(n, min, max, values, prob=rep(1,length(values)))
#ISALIAS dempiricalC
#--------------------------------------------
{ 
  if(length(n) > 1) n <- length(n)
  return(qempiricalC(runif(n),min=min,max=max,values=values,prob=prob))

  }
