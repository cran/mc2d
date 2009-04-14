#<<BEGIN>>
Ops.mcnode <- function(e1,e2)
#TITLE Operations on mcnode Objects
#DESCRIPTION
# This function alters the way operations are performed on \samp{mcnode} objects for a better consistancy of the theory.
#KEYWORDS utilities
#INPUTS
#{e1}<<An \samp{mcnode} object, a vector or an array.>>
#{e2}<<An optionnal \samp{mcnode} object, a vector or a matrix with at least one of both objects as an \samp{mcnode}.>>
#VALUE
#The results as a \samp{mcnode} object.
#DETAILS
#This method will be used for any of the Group \code{\link{Ops}} functions.
#
#The rules are as following (illustrated with a \samp{+} function and ignoring the \samp{nvariates} dimension):
#{*}<<\samp{0 + 0 = 0};>>
#{*}<<\samp{0 + V = V}: classical recycling of the scalar;>>
#{*}<<\samp{0 + U = U}: classical recycling of the scalar;>>
#{*}<<\samp{0 + VU = VU}: classical recycling of the scalar;>>
#{*}<<\samp{V + V = V}: if both of the same \samp{(nsv)} dimension;>>
#{*}<<\samp{V + U = VU}: the \samp{U} object will be recycled "by row". The \samp{V} object will be recycled classically "by column";>>
#{*}<<\samp{V + VU = VU}: if the dimension of the \samp{V} is \samp{(nsv)} and the dimension of the \samp{VU} is \samp{(nsv x nsu)}. The \samp{V} object will be recycled classically "by column";>>
#{*}<<\samp{U + U = U}: if both of the same \samp{(nsu)} dimension;>>
#{*}<<\samp{U + VU = VU}: if the dimension of the \samp{U} is \samp{(nsu)} and the dimension of the \samp{VU} is \samp{(nsv x nsu)}. The \samp{U} object will be recycled "by row";>>
#{*}<<\samp{VU + VU = VU}: if the dimension of the \samp{VU} nodes is \samp{(nsu x nsv)};>>
#
#A vector or an array may be combined with an \samp{mcnode} of size \samp{(nsv x nsu)} if an \samp{mcnode} of this dimension
#may be built from this vector/array using the \samp{mcdata} function. See \code{\link{mcdata}} for the rules.
#
#The \samp{outm} attribute is transferred as following: \samp{each + each = each}; \samp{none + other = other};
#\samp{other1 + other2 = other1}. The \samp{outm} attribute of the resulting node may be changed using the \code{\link{outm}} function.
#
#For multivariate nodes, a recycling on the \samp{nvariates} dimension is done if a \samp{(nsu x nsv x nvariates)} node
#is combined with a \samp{(nsu x nsv x 1)} node.

#SEE ALSO
#\code{\link{mcdata}}, \code{\link{mcstoc}}
#EXAMPLE
#oldvar <- ndvar()
#oldunc <- ndunc()
#ndvar(30)
#ndunc(20)
#
### Given
#x0 <- mcdata(3,type="0")
#xV <- mcdata(1:ndvar(),type="V")
#xU <- mcdata(1:ndunc(),type="U")
#xVU <- mcdata(1:(ndunc()*ndvar()),type="VU")
#x0M <- mcdata(c(5,10),type="0",nvariates=2)
#xVM <- mcdata(1:(2*ndvar()),type="V",nvariates=2)
#xUM <- mcdata(1:(2*ndunc()),type="U",nvariates=2)
#xVUM <- mcdata(1:(2*(ndunc()*ndvar())),type="VU",nvariates=2)
#
### All possible combinations
### "0"
#-x0
#x0 + 3
#
### "V"
#-xV
#3 + xV
#xV * (1:ndvar())
#xV * x0
#xV - xV
#
### "U"
#-xU
#xU + 3
#(1:ndunc()) * xU
#xU * x0
#xU - xU
#
### Watch out the resulting type
#xV + xU
#xU + xV
#
### "VU"
#-xVU
#3 + xVU
#(1:(ndunc()*ndvar())) * xVU
#xVU + xV
#x0 + xVU
#xU + xVU
#xVU - xVU
#
### Some Multivariates
#x0M+3
#xVM * (1:ndvar())
#xVM - xV
#xUM - xU
#xVUM - xU
#AUTHOR Regis Pouillot
#CREATED 08-01-25
#--------------------------------------------
{
  err <- "Incompatible mcnode dimensions"

  if(missing(e2)) {
    type <- attr(e1,"type")
    outm <- attr(e1,"outm")}   # only e1 like -X

  else {
      if(!inherits(e1,"mcnode")) {   # only e2 mcnode
        dime <- dim(e2)
        type <- attr(e2,"type")
        e1 <- mcdata(e1,type=type, nsv=dime[1], nsu=dime[2],nvariates=dime[3],outm="none")
        }
      else if(!inherits(e2,"mcnode")) {   # only e2 mcnode
        dime <- dim(e1)
        type <- attr(e1,"type")
        e2 <- mcdata(e2,type=type, nsv=dime[1], nsu=dime[2],nvariates=dime[3],outm="none")
        }

      dim1 <- dim(e1)
      dim2 <- dim(e2)
      dimf <- pmax(dim1,dim2)

      outm1 <- attr(e1,"outm")
      outm2 <- attr(e2,"outm")
      if(is.null(outm1) && is.null(outm2)) outm <- "each"
      else if(is.null(outm1)) outm <- outm2
      else if(is.null(outm2)) outm <- outm1
      else if(outm1 == "each" || outm2 == "each") outm <- "each"
        else if(outm1=="none") outm <- outm2
          else outm <- outm1

      type1 <- attr(e1,"type")
      type2 <- attr(e2,"type")
      if(type1==type2){ type <- type1 ; if(any(dim1[1:2] != dim2[1:2])) stop(err)}        # U+U, V+V, VU + VU, gère les deux premieres dimensions

      else if(type1=="0"){ type <- type2 ; e1 <- array(rep(e1,each=dimf[1]*dimf[2]),dim=dimf)}
      else if(type2=="0"){ type <- type1 ; e2 <- array(rep(e2,each=dimf[1]*dimf[2]),dim=dimf)}

      else {                                                                              #U+V, V+U, V+VU and U+VU
        type <- "VU"

        if(type1=="U") {                                                                    #U+V and U+VU
          if(type2=="VU" && dim1[2]!= dim2[2]) stop(err)
          e1 <- array(apply(e1,3,matrix,ncol=dimf[2],nrow=dimf[1],byrow=TRUE),dim=dimf)     # recycling dim 1,2,3
          dim1 <- dimf}
        else if(type2=="U") {                                                               #V+U and VU+U dim 1,2,3
          if(type1=="VU" && dim1[2]!= dim2[2]) stop(err)
          e2 <- array(apply(e2,3,matrix,ncol=dimf[2],nrow=dimf[1],byrow=TRUE),dim=dimf)     # recycling dim 1,2,3
          dim2 <- dimf}
        if(dim1[1]!= dim2[1]) stop(err)                                                     # reste V+VU, VU+V, no problem of recycling
        if(type1=="V") e1 <- array(apply(e1,3,matrix,ncol=dimf[2],nrow=dimf[1]),dim=dimf)  # necessary recycling of V for arrays
        else if(type2 == "V") e2 <- array(apply(e2,3,matrix,ncol=dimf[2],nrow=dimf[1]),dim=dimf)
      }

        if(dim1[3] != dimf[3]){                      # gère la troisieme dimension
          if(dim1[3] != 1) stop(err)
          else e1 <- array(e1,dim=dimf)}

        else if(dim2[3] != dimf[3]){
          if(dim2[3]!=1) stop(err)
          else e2 <- array(e2,dim=dimf)}
    }

    res <- NextMethod(.Generic)
    class(res) <- "mcnode"
    attr(res,"type") <- type
    attr(res,"outm") <- outm
    return(res)
}

