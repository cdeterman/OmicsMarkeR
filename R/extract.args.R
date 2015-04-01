## Function to extract arguments from previously fit model
#' @title Argument extractor
#' @description Extract arguments from previously fs.stability models
#' @param fs.model Previously fit fs.stability model
#' @param method Which model to extract from
#' @return args List of model arguments

extract.args <- function(fs.model, method)
{
    ind <- which(names(fs.model$performance) == method)
    metrics <- data.frame(fs.model$performance[[ind]])
    
    args <- switch(method,
                   plsda = {out <- list(ncomp = metrics$ncomp)},
                   gbm = {out <- c(metrics[,colnames(metrics) %in% 
                                               c("interaction.depth", 
                                                 "n.trees", "shrinkage")])},
                   rf = {out <- list(mtry = metrics$mtry)},
                   svm = {out <- list(C = metrics$C)},
                   pam = {out <- list(threshold = metrics$threshold)},
                   glmnet = {out <- c(metrics[colnames(metrics) %in% 
                                                  c("alpha", "lambda")])}
    )
    return(args)
}