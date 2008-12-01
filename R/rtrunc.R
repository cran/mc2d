#<<BEGIN>>
rtrunc <- function(distr=runif, n, linf=-Inf, lsup=Inf,...)
#TITLE Random Truncated Distributions
#DESCRIPTION
# Provides samples from classical \R distributions and \samp{mc2d} specific
# distributions truncated between \samp{linf} and \samp{lsup}.
#KEYWORDS distribution
#INPUTS
#{distr}<<A function providing random data or its name as character.
#The function 'rdistr' should have a 'qdistr' form (with argument 'p') and a 'pdistr' form
#(with argument 'q'). Example : 'rnorm' (has a 'qnorm' and a 'pnorm' form), 'rbeta', 'rbinom', 'rgamma', ...>>
#{n}<<The size of the sample.>>         .
#[INPUTS]
#{linf}<<A vector of lower bounds.>>
#{lsup}<<A vector of upper bounds.>>
#{\dots}<<All arguments to be passed to \samp{pdistr} and \samp{qdistr}.>>
#VALUE
#A vector of \samp{n} values.
#DETAILS
#The function 1) evaluates the \samp{p} values corresponding to \samp{linf} and \samp{lsup} using \samp{pdistr};
#2) samples \samp{n} values using \samp{runif(n, min=pinf, max=psup)}, and 3) takes
#the \samp{n} corresponding quantiles from the specified distribution using \samp{qdistr}.
#
#All distributions (but sample) implemented in the stats library could be used.
#The arguments in \dots should be named. Do not use 'log' or 'log.p' or 'lower.tail'.
#NOTE
#The inversion of the quantile function leads to time consuming functions for some distributions.
#EXAMPLE
#rtrunc("rnorm", n=10, linf=0)
#range(rtrunc(rnorm, n=1000, linf=3, lsup=5, sd=10))
#AUTHOR Regis Pouillot
#CREATED 08-02-20
#--------------------------------------------
{
    if(!is.character(distr)) distr <- as.character(match.call()$distr)          #retrieve the name of the function
    distr <- substr(distr, 2, 1000)                                             #remove the r

    if(any(pmin(linf,lsup)!=linf)) stop("linf should be <= lsup")  #recycle vectors

    pfun <- get(paste("p",distr,sep=""),mode="function")

    pinf <- pfun(q=linf,...)
    psup <- pfun(q=lsup,...)

    p <- runif(n,min=pinf,max=psup)

    qfun <- get(paste("q",distr,sep=""),mode="function")

    res <- qfun(p,...)
    # Some possible problem (check if you think to others)
    #
    res[pinf==0 & res > lsup] <- NaN          #ex: rtrunc("lnorm",10,linf=-2,lsup=-1)
    res[psup==1 & res < linf] <- NaN          #ex: rtrunc("unif",10,linf=2,lsup=4,max=1)
    res[is.na(linf) | is.na(lsup)] <- NaN   #ex: rtrunc("norm",10,sd=-2)

    return(res)
}
