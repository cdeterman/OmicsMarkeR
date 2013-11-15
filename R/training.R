
#' @title Model Training
#' @description This fits each model with the defined parameters
#' @param data Dataframe consisting of both numeric feature values and a single column
#' named '.classes' to denoted group membership.
#' @param method String dictating which model to fit
#' @param tuneValue List of parameters to be applied to the specific model
#' @param obsLevels Observed group levels
#' @param theDots List of additional parameters to be applied to the specific model
#' @return \item{fit}{Fitted model with list with the following elements:}
#' @return \itemize{
#'  \item{xNames: Names of the features}
#'  \item{tuneValue: Parameters applied to the fitted model}
#'  \item{obsLevels: Observed levels of the groups}}
#' @author Charles Determan Jr
# @import DiscriMiner
#' @import randomForest
#' @import caret
#' @import e1071
#' @import gbm
#' @import pamr
#' @import glmnet
#' @export

training <-
  function(data, method, tuneValue, obsLevels, theDots = NULL)
    {
    
    data <- as.data.frame(data)

    ## pam and will crash if there is a resample with <2 observations
    ## in a class. We will detect this and remove those classes.
    if(method == "pam")
    {
      yDist <- table(data$.classes)
      if(any(yDist < 2))
      {
        smallClasses <- names(yDist[yDist < 2])
        data <- data[!(data$.classes %in% smallClasses),]
      }
    }
    
    ## Factor the class labels
    data$.classes <- factor(as.character(data$.classes), levels = obsLevels)
    xNames <- names(data)[!(names(data) %in% ".classes")]
    
    trainX <- data[,!(names(data) %in% ".classes"), drop = FALSE]
    trainY <- data[,".classes"]
    
    if(method == "gbm" & length(obsLevels) == 2)  numClasses <- ifelse(data$.classes == obsLevels[1], 1, 0)
                       
    modelFit <- switch(method,
                       plsda =
                         {
                           #library(DiscriMiner)
                           # retain.models omitted because when this is used, the final model is only using the best component
                           # may switch to retain all models but this omits more processing that is likely superfluous
                           
                           # check for number of components provided.  This is important following selection of the best model
                           if(tuneValue$.ncomp == 1){
                             warning("PLSDA model contained only 1 component. PLSDA requires at least 2 components.\nModel fit with 2 components")
                             tuneValue$.ncomp = 2
                           }
                           
                           plsDA(trainX, 
                                 trainY,
                                 autosel=F,
                                 validation = NULL,
                                 comps = tuneValue$.ncomp,
                                 cv ="none")
                         },
                       
                       gbm =  
                         {
                           #library(gbm)
                           # need to make sure only extract arguments that pertain to gbm
                           gbm.args <- c("w", "var.monotone", "n.minobsinnode", 
                                         "bag.fraction", "var.names", "response.name", "group") 
                           theDots <- theDots[names(theDots) %in% gbm.args]
                           
                           if(ncol(trainX) < 50 | nrow(trainX) < 50){
                             if(is.null(theDots) | length(theDots) == 0){
                               if(nrow(trainX) < 30){
                                 theDots <- list(n.minobsinnode = 2)
                               }else{
                                 theDots <- list(n.minobsinnode = 5)  
                               }
                             }
                           }
                           
                           # determine if binary or multiclass
                           gbmdist <- if(length(unique(trainY)) == 2){
                             "bernoulli"}else{
                               "multinomial"
                             }         
                           
                           # check gbm setup file to see if this is necessary
                           modY <- if(gbmdist != "multinomial") numClasses else trainY
                           
                           if(gbmdist != "multinomial"){
                             modY <- numClasses
                           }else{
                             modY <- trainY
                           }
                           
                           modArgs <- list(x = trainX,
                                           y = modY,
                                           interaction.depth = tuneValue$.interaction.depth,
                                           n.trees = tuneValue$.n.trees,
                                           shrinkage = tuneValue$.shrinkage, 
                                           distribution = gbmdist,
                                           verbose = FALSE)
                           
                           if(length(theDots) > 0) modArgs <- c(modArgs, theDots)
                           
                           do.call("gbm.fit", modArgs)
                           
                         },
                       
                       rf =
                         {
                           #library(randomForest)
                           
                           rf.args <- c("maxnodes", "keep.forest", "keep.inbag")
                           theDots <- theDots[names(theDots) %in% rf.args]
                           
                           modArgs <- list(x = trainX,
                                           y = trainY,
                                           importance = TRUE,
                                           mtry = tuneValue$.mtry,
                                           ntree=round.multiple(sqrt(ncol(trainX)), target = 50)
                           )
                           
                           if(length(theDots) > 0) modArgs <- c(modArgs, theDots)
                           
                           do.call("randomForest", modArgs)
                         },
                       
                       svm =
                         { 
                           #library(e1071)
                           out <- svm(trainX,
                                      trainY,
                                      cost = tuneValue$.C, 
                                      cachesize=500,
                                      scale=F, 
                                      type="C-classification", 
                                      kernel="linear")                         
                           out
                         },
                       
                       pam = 
                         {
                           #library(pamr)    
                           pamr.args <- c("n.threshold", "threshold.scale", "scale.sd", "se.scale")
                           theDots <- theDots[names(theDots) %in% pamr.args]
                           
                           modArgs <- list(data = list(x = t(trainX), y = trainY, geneid = as.character(colnames(trainX))),
                                           threshold = tuneValue$.threshold
                           )
                           
                           if(length(theDots) > 0) modArgs <- c(modArgs, theDots)
                           
                           do.call("pamr.train", modArgs)
                         }, 
                       
                       glmnet =
                         {
                           #library(glmnet)
                           numLev <- if(is.character(trainY) | is.factor(trainY)) length(levels(trainY)) else NA
                           
                           glmnet.args <- c("offset", "nlambda", "weights", "standardize","intecept", 
                                            "dfmax", "pmax","exclude","penalty.factor","lower.limits",
                                            "upper.limits","maxit","standardize.response","type.multinomial")
                           theDots <- theDots[names(theDots) %in% glmnet.args]
                           
                           if(!is.null(theDots)){
                             if(all(names(theDots) != "family"))
                             {
                               if(!is.na(numLev))
                               {
                                 fam <- ifelse(numLev > 2, "multinomial", "binomial")
                               } else stop("Error: levels of classes couldn't be determined for glmnet")
                               
                               if(is.null(theDots)){
                                 theDots <- list(family = fam)
                               }else{
                                 theDots$family <- fam
                               }
                             }
                           }else{
                             if(!is.na(numLev))
                             {
                               fam <- ifelse(numLev > 2, "multinomial", "binomial")
                             } else stop("Error: levels of classes couldn't be determined for glmnet")
                             
                             if(is.null(theDots)){
                               theDots <- list(family = fam)
                             }
                           }
                           
                           modelArgs <- c(
                             list(x = as.matrix(trainX), y = trainY, alpha = tuneValue$.alpha),
                             theDots)
                           
                           out <- do.call("glmnet", modelArgs) 
                           out 
                         }
                     )
    
    modelFit$xNames <- xNames
    modelFit$tuneValue <- tuneValue
    modelFit$obsLevels <- obsLevels
    
    list(fit = modelFit)
}
