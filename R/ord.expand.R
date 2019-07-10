

# ------------------------------------------------------------
# function to expand ordinal scores and explanatory variates
# ------------------------------------------------------------

ord.expand <- function(scores,data,subjects,categories){
  
  #attach(data,warn.conflicts=FALSE)
  namvars <- names(data)
  nvars <- length(names(data))
  categs1 <- categories-1
  exdsize <- categs1*dim(data)[1]
  exdata <- matrix(0,ncol=nvars,nrow=exdsize)
  for (j in 1:nvars){
    if (namvars[j]==scores){
      boscores <- matrix(0,ncol=categs1,nrow=dim(data)[1])
      for (i in 1:categs1){
        boscores[,i] <- as.numeric(data[,j]<=i)
      } # for i
      exdata[,j] <- as.vector(t(boscores))
      cuts <- factor(rep(1:categs1,dim(data)[1]))
    }
    else{
      exdata[,j] <- rep(data[,j],each=categs1)
    } # if
    if (namvars[j]==subjects) namvars[j] <- "id"
  } # for j
  exdata <- as.data.frame(exdata,row.names=1:exdsize)
  names(exdata) <- namvars
  exdata <- cbind(cuts,exdata)
  #detach(data)
  
  # output
  list(exdata=exdata)
}