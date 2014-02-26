
#' @title Class Prediction
#' @description This function evaluates a single fitted model and returns the predicted group memberships of new data.
#' @param method String of the model to be evaluated
#' @param modelFit The fitted model being evaluated
#' @param orig.data The orginal data before subsetting training sets.  Required to have the 'observed' group membership
#' @param indicies The indicies for the training subsets
#' @param newdata The testing data to predict group membership
#' @param parms Optional alternate parameters being fit to the model
#' @return Returns a list of predicted group membership
#' @import DiscriMiner
#' @import randomForest
#' @import caret
#' @import e1071
#' @import gbm
#' @import pamr
#' @import glmnet
#' @import data.table
#' @export

predictNewClasses <- function(modelFit, method, orig.data, newdata, parms = NULL)
{
  # check for which column contains classes
  o.factors <- as.vector(which(sapply(orig.data, is.factor)))
  if(length(o.factors) > 1){
    stop("\n Error: your data contains multiple factor columns.")
  }else{
    o.name <- colnames(orig.data[,o.factors])
    colnames(orig.data)[o.factors] = ".classes"
    obsLevels <- levels(orig.data[,o.factors])
  }
  
  # check for any factors in new data
  if(any(sapply(newdata,is.factor))){
    warning("\n Columns identified with factors have been removed from newdata.")
    f.rm <- as.vector(which(sapply(newdata, is.factor)))
    newdata[,f.rm] <- NULL}
  
  # extract model to be evaluated
  focusModel <- which(names(modelFit$performance) == method)
  modelFit <- modelFit$performance[[focusModel]]
  
  # extract model parameters
  pars <- sapply(params(method), FUN = function(x) x$parameter)
  tVal <- modelFit[,c(as.vector(pars))]
  names(tVal) <- paste(".", as.vector(pars), sep = "")
  tVal <- as.data.frame(t(tVal))

  # get trained model
  train.model <- training(data = orig.data, method, tuneValue = tVal, obsLevels = obsLevels, theDots = parms)

  predictedValue <- switch(method,                           
                           plsda =
                             {
                               # library(DiscriMiner)
                               # check for number of components provided.  This is important following selection of the best model
                               ncomp <- tVal
                               if(ncomp == 1){
                                 warning("PLSDA model contained only 1 component. PLSDA requires at least 2 components.\nModel fit with 2 components")
                                 ncomp = 2
                               }                        
                               
                               o.nr <- seq(nrow(orig.data))
                               n.nr <- seq(from = (nrow(orig.data)+1), to = nrow(newdata) + nrow(orig.data))
                               full.data <- as.data.frame(rbindlist(list(orig.data, newdata)))
                               
                               tmp <- plsDA(full.data[,-which(names(full.data) %in% c(".classes"))], 
                                            full.data[,c(".classes")],
                                            autosel=F,
                                            learn = o.nr,
                                            test = n.nr,
                                            validation = "learntest",
                                            comps = ncomp,
                                            cv ="none",
                                            retain.models = TRUE)$classification
                               
                               if(tVal < 2){
                                 out <- lapply(tmp, as.character)[[1]]
                               }else{
                                 out <- lapply(tmp, as.character)
                               }
                               
                               out
                             },

                           gbm =
                             {
                               if(is.null(parms)){
                                 #library(gbm)                               
                                 gbmProb <- predict(train.model$fit, 
                                                    newdata=newdata[,!names(newdata) == ".classes"], 
                                                    type = "response",
                                                    n.trees = train.model$fit$n.trees)
                                 gbmProb[is.nan(gbmProb)] <- NA
                                 
                                 # need a check if all NA
                                 # if so, n.trees are way too high                               
                                 if(length(obsLevels) <= 2)
                                 {
                                   out <- ifelse(gbmProb >= .5, obsLevels[1], obsLevels[2])
                                   ## to correspond to gbmClasses definition above
                                 } else {
                                   out <- colnames(gbmProb)[apply(gbmProb, 1, which.max)]
                                 }
                               }else{
                                 tmp <- predict(train.model$fit, 
                                                newdata[,!names(newdata) == ".classes"], 
                                                type = "response", 
                                                n.trees = parms$n.trees)
                                 
                                 if(length(obsLevels) <= 2){
                                   # if only one other parameter, need to convert to matrix
                                   if(is.vector(tmp)) tmp <- matrix(tmp, ncol = 1)
                                   tmp <- apply(tmp, 2,
                                                function(x, nm = obsLevels) ifelse(x >= .5, nm[1], nm[2]))
                                 }else{
                                   tmp <- apply(tmp, 3,
                                                function(y, nm = obsLevels) nm[apply(y, 1, which.max)])
                                 }
                                 
                                 # convert to list compatible splits
                                 if(!is.list(tmp)) tmp <- split(tmp, rep(1:ncol(tmp), each = nrow(tmp)))
                                 out <- as.vector(unlist(tmp))
                               }
                               out
                             },
                           
                           rf =
                             {
                              #library(randomForest)
                               out <-  as.character(predict(train.model, newdata))
                               out
                             },
                           svm =                           
                             {
                               #library(e1071)                             
                               out <- as.character(predict(train.model, newdata = newdata))
                               out
                             },
                           pam =
                             {
                               #library(pamr)
                               out <- as.character(
                                 pamr.predict(train.model,
                                              t(newdata),
                                              threshold = train.model$tuneValue$.threshold))
                               #pamr.predict
                               if(!is.null(param))
                               {
                                 tmp <- vector(mode = "list", length = nrow(param) + 1)
                                 tmp[[1]] <- out
                                 for(j in seq(along = param$.threshold))
                                 {
                                   tmp[[j+1]] <- as.character(
                                     pamr.predict(
                                       train.model,
                                       t(newdata),
                                       threshold = param$.threshold[j]))
                                 }
                                 out <- tmp
                               }
                               out
                             },
                           
                           glmnet =
                             {                          
                               #library(glmnet)
                               if(!is.matrix(newdata)) newdata <- as.matrix(newdata)
                               
                               if(!is.null(param))
                               {
                                 out <- predict(train.model, newdata, s = param$.lambda, type = "class")
                                 out <- as.list(as.data.frame(out, stringsAsFactors = FALSE))
                               } else {
                                 
                                 if(is.null(train.model$lambdaOpt))
                                   stop("optimal lambda not saved; needs a single lambda value")
                                 
                                 out <- predict(train.model, newdata, s = train.model$lambdaOpt, type = "class")[,1]
                               }
                               out
                             },
                           
  )
  data.frame(predictedClass = predictedValue)
}

