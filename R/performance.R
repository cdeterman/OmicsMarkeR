
#' @title Performance Statistics (Internal for \code{perf.calc})
#' @description Calculates confusion matrix and ROC statistics comparing the results of the fitted models
#' to the observed groups.
#' @param pred vector of groups predicted by a fitted classification model
#' @param obs vector of groups from the original dataset
#' @return Returns confusion matrix and ROC performance statistics including
#' Accuracy, Kappa, ROC.AUC, Sensitivity, Specificity, Positive Predictive Value, and Negative Predictive Value
#' @seealso caret function \code{\link{confusionMatrix}}
#' @import caTools
#' @importFrom caret confusionMatrix


performance.stats <- function(pred, obs)
{
  isNA <- is.na(pred)
  pred <- pred[!isNA]
  obs <- obs[!isNA]
  
  if(length(obs) + length(pred) == 0)
  {
    out <- rep(NA, 2)
  } else {
    #library(caTools)
    tmp.auc <- colMeans(colAUC(order(pred), obs, plotROC = FALSE, alg = "ROC"))
    pred <- factor(pred, levels = levels(obs))  
    
    tmp <- confusionMatrix(pred, obs)
    
    overall.stats <- c("Accuracy", "Kappa")
    byclass.stats <- c("Sensitivity", "Specificity", "Pos Pred Value", "Neg Pred Value")
    
    colnames(tmp$byClass)
    
    overall.index <- which(names(tmp$overall) %in% overall.stats)
    if(nlevels(obs) == 2){
      byclass.index <- which(names(tmp$byClass) %in% byclass.stats)
    }else{
      byclass.index <- which(colnames(tmp$byClass) %in% byclass.stats)
    }
    
    out <- c(tmp$overall[overall.index], dim(tmp$table)[1], tmp$byClass[byclass.index], tmp.auc)
  }
  names(out) <- c("Accuracy", "Kappa", "Table", "Sensitivity", "Specificity", "Pos Pred Value", "Neg Pred Value", "ROC.AUC")
  out <- out[c("Accuracy", "Kappa", "ROC.AUC", "Sensitivity", "Specificity", "Pos Pred Value", "Neg Pred Value", "Table")]
  
  if(any(is.nan(out))) out[is.nan(out)] <- NA
  out
}

#' @title Performance Statistics Calculations
#' @description Calculates confusion matrix and ROC statistics comparing the results of the fitted models
#' to the observed groups.
#' @param data dataframe of predicted (pred) and observed (obs) groups
#' @param lev Group levels
#' @param model String indicating which model was initially run
#' @return Returns confusion matrix and ROC performance statistics including
#' Accuracy, Kappa, ROC.AUC, Sensitivity, Specificity, Positive Predictive Value, and Negative Predictive Value
#' @seealso caret function \code{\link{confusionMatrix}}
# ' @export

perf.calc <- function(data,          
                      lev = NULL, 
                      model = NULL)     
{
  if(is.character(data$obs)){
    data$obs <- factor(data$obs, levels = lev)
    data$pred <- factor(data$pred, levels = lev)
  } 
  out <- performance.stats(data[,"pred"], data[,"obs"])
  out <- out[!names(out) == "Table"]
  out
}
