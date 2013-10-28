

#' @title Performance Metrics of fs.stability or fs.ensembl.stability object
#' @description This will provide a concise data.frame of confusion matrix and ROC statistics
#' from the results of \code{fs.stability} or \code{fs.ensembl.stability}.
#' @param fit.model An fs.stability or fs.ensembl.stability object
#' @param digits How many digits to round values
#' @return Dataframe of performance statistics by model
#' @author Charles E. Determan Jr.
#' @export

performance.metrics <- function(fit.model,   # fs.stability object
                                digits = max(3, getOption("digits") - 3))
{
  cat("\nModel Performance Statistics\n\n") 
  
  # need to set class of performance statistics in fs.stability object
  # instead, probably be easiest to simply use the fs.stability results and draw performance with generic performance.metrics function
 
  perf <- fit.model$performance
  fit.models <- fit.model$methods
  
  if(length(fit.models) >= 2){
    mult.param <- sapply(fit.models, function(x) as.character(params(x)[[1]]$parameter))
    out <- mapply(perf, FUN = function(x,y) x[,!colnames(x) %in% y], y = mult.param)
    out
  }else{
    perf <- do.call(rbind, perf)
    param <- as.character(params(fit.models)[[1]]$parameter)
    tmp <- perf[,!colnames(perf) %in% param]
    out <- as.data.frame(tmp)
    colnames(out) <- fit.models
    out
  }
}