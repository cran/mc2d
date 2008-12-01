#<<BEGIN>>
plot.mc <- function(x, prec=0.01, stat = c("median","mean"), lim = c(0.025,0.975), na.rm=TRUE, griddim = NULL, xlab = NULL, ylab = "Fn(x)", main = "", draw = TRUE, ...)
#TITLE Plots Results of a Monte Carlo Simulation
#DESCRIPTION
# Plots the empirical cumulative distribution function of a \samp{mcnode} or a \samp{mc} object ("0" and "V" nodes) or the
#empirical cumulative distribution function of the estimate of a \samp{mcnode} or \samp{mc} object ("U" and "VU" nodes).
#KEYWORDS hplot
#INPUTS
# {x}<<a \samp{mcnode} or a \samp{mc} objects>>
#[INPUTS]
#{prec}<<the precision of the plot. 0.01 will
#provide an ecdf from the 0.00, 0.01, .02,  ..., 1.00 quantiles, 0.001 will provide
#a 0.000, 0.001, 0.002, ..., 1.000 quantiles,... >>
#{stat}<<the function used for estimates (2D \samp{mc} or \samp{mcnode}). By default the median.>>
#{lim}<<a vector of numbers (between 0 and 1) indicating the enveloppe (2D \samp{mc} or \samp{mcnode}) . Maybe \samp{NULL} or empty.>>
#{na.rm}<<Should NA values be discarded>>
#{griddim}<<a vector of two integers, indicating the size of the grid of the graph. If \samp{NULL}, the grid is calculated to produce a "nice" graph.>>
#{xlab}<<vector of labels for the x-axis. If \samp{NULL}, use the name of the node.>>
#{ylab}<<vector of labels for the y-axis.>>
#{main}<<vector of main titles of the graph.>>
#{draw}<<Should the plot be drawn?>>
#{\dots}<<further arguments to be passed to \samp{plot.ecdf}.>>
#DETAILS
#\samp{plot.mcnode} is a user-friendly function that send the \samp{mcnode} to \samp{plot.mc}.</>
#For \samp{"VU"} and \samp{"U"} \samp{mcnode}s, quantiles are calculated using \code{\link{quantile.mc}}
#within each of the \samp{nsu} simulations (i.e. by columns of each \samp{mcnode}). The medians (but may be
#the means using \samp{stat="mean"}) calculated from the \samp{nsu} values are plotted. The 0.025 and 0.975 quantiles (default values
#of \samp{lim}) of these quantiles are used as the enveloppe.
#REFERENCE
#Cullen AC and Frey HC (1999) Probabilistic techniques in exposure assessment. Plenum Press, USA, pp. 81-155. 
#VALUE
#A \samp{plot.mc} object, list of the quantiles used to plot the draw.
#SEE ALSO
#\code{\link{ecdf}}, \code{\link{plot}}, \code{\link{quantile.mc}}
#EXAMPLE
#data(total)
#plot(xVUM3)
#plot(total)
#AUTHOR Regis Pouillot
#CREATED 07-08-01
#REVISED 07-08-01
#--------------------------------------------
#
{
  if(inherits(x,"mc")==TRUE) {
    x <- quantile.mc(x, probs=seq(0,1,prec),lim = lim, na.rm=na.rm, lnames=xlab)
    }

  if(draw) {
    y <- x                           # for a correct return
    stat <- match.arg(stat)

	 beau <- function(n){
		nc <- round(sqrt(n))
		nr <- ceiling(n/nc)
		c(nc,nr)}

   noms <- names(rapply(x,function(x) 1))    #moche mais efficace
   if(is.null(xlab)) xlab <- noms
   n <- length(noms)

	 main <- rep(main,n)
	 xlab <- rep(xlab,n)
	 ylab <- rep(ylab,n)

  if(is.null(griddim)) griddim <- beau(n)
  if(prod(griddim) < n) op <- par(mfrow=griddim,ask=TRUE)
     else op <- par(mfrow=griddim )

  try({   #to restore par in case of error

  i <- 1
  env <- environment()
  
  LEPLOT <- function(y,...){
      if(nrow(y) != 1) {
        if(stat=="median") y <- y[-2,,drop=FALSE]
        else y <- y[-1,,drop=FALSE]}                                              #Retrait median or mean
  		nr <- nrow(y)
      xlima <- range(y,na.rm=TRUE)
  		if(xlima[1]==xlima[2]) xlima <- xlima + c(-0.5,0.5)
      i <- get("i",envir=env)
      plot.ecdf(y[1,],main=main[i],xlim=xlima,ylab = ylab[i], verticals = TRUE, do.points = FALSE, xlab=xlab[i],...)
      if(nr > 1){
        for(j in 2:nr) plot.ecdf(y[j,],verticals=TRUE, do.points=FALSE, add=  TRUE, lty=2,...)
      }
    assign("i",i+1,envir=env) }
    
  rapply(y,LEPLOT)
    
  })
	par(op)
  }
  class(x) <- "plotmc"
  return(invisible(x))}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<BEGIN>>
plot.mcnode <- function(x, ...)
#ISALIAS plot.mc
#--------------------------------------------
{ nom <- deparse(substitute(x))
  x <- list(x)
  names(x) <- nom
  class(x) <- "mc"
  x <- plot.mc(x, ... )
  return(invisible(x))}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<BEGIN>>
plot.plotmc <- function(x, ...)
#ISALIAS plot.mc
#--------------------------------------------
{ x <- plot.mc(x, ... )
  return(invisible(x))}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

