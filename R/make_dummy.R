# -----------------------------------------------------------------
# function to create dummy variables for outcome and predictors
# -----------------------------------------------------------------

#' @export

make.dummy <- function(data, resp.ind, nresponse, xvars.expand=NULL, xvars.common=NULL){
  
  data.dummy <- eval(data)
  resp.dummy <- dummy_cols(data.dummy$resp.ind)
  names(resp.dummy) <- c("resp.ind", paste0("resp.ind", 1:nresponse))
  resp.dummy1 <- resp.dummy[,-c(1:2)]
  
  if(is.null(xvars.expand)){
    
    commonxvars <- as.data.frame(data.dummy[,xvars.common])
    names(commonxvars) <- xvars.common
    cbind(resp.dummy, commonxvars)
    
  }else{

    nexpand <- length(xvars.expand)
    
    xvars.dummy.list <- list()
    for (i in 1:nexpand){
      xvars.dummy <- resp.dummy1 * data.dummy[,xvars.expand[i]]
      xvars.dummy <- cbind(data.dummy[,xvars.expand[i]], xvars.dummy)
      names(xvars.dummy) <- c(xvars.expand[i], paste0(xvars.expand[i], 2:nresponse))
      xvars.dummy.list[[i]] <- xvars.dummy
    }
    
    xvars.dummy.frame <- do.call("cbind", xvars.dummy.list)
    
    if (is.null(xvars.common)){
      cbind(resp.dummy, xvars.dummy.frame)
    }else{
      xvars.common.frame <- data.frame(data.dummy[,xvars.common])
      names(xvars.common.frame) <- xvars.common
      cbind(resp.dummy, xvars.dummy.frame, xvars.common.frame)
    }
  }

}

#dummyout <- make.dummy(data=out[[3]], resp.ind, nresponse=out[[2]], xvars.expand=out[[8]], xvars.common=out[[9]])
#dummyout
