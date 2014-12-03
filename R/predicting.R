
#' @title Model Group Prediction
#' @description This function evaluates a single fitted model and returns the predicted group memberships.
#' @param method String of the model to be evaluated
#' @param modelFit The fitted model being evaluated
#' @param orig.data The orginal data before subsetting training sets.  Required to have the 'observed' group membership
#' @param indicies The indicies for the training subsets
#' @param newdata The testing data to predict group membership
#' @param param The parameters being fit to the model (Determined by model optimization).
#' @return Returns a list of predicted group membership
#' @import DiscriMiner
#' @import randomForest
#' @import e1071
#' @import gbm
#' @import pamr
#' @import glmnet
# ' @export


predicting <- function(method, modelFit, orig.data, indicies, newdata, param = NULL)
{
  if(any(colnames(newdata) == ".classes")) newdata$.classes <- NULL

  coerceChar <- function(x){
    as.data.frame(lapply(x, as.character), stringsAsFactors = FALSE)
  }  
  
  predictedValue <- switch(method,                           
                           plsda =
                           {
                             # require(DiscriMiner)
                             # check for number of components provided.  This is important following selection of the best model
                             
                             ncomp <- modelFit$tuneValue$.ncomp
                             if(ncomp == 1){
                               warning("PLSDA model contained only 1 component. PLSDA requires at least 2 components.\nModel fit with 2 components")
                               ncomp = 2
                             }                        
                             
                             vars <- as.matrix(orig.data[,-which(names(orig.data) %in% c(".classes"))])
                             mode(vars) <- 'numeric'
                             
                             tmp <- plsDA(vars, 
                                          orig.data[,c(".classes")],
                                          autosel=F,
                                          learn = indicies,
                                          test = seq(nrow(orig.data))[-unique(indicies)],
                                          validation = "learntest",
                                          comps = ncomp,
                                          cv ="none",
                                          retain.models = TRUE)$classification
                             
                             if(ncomp < 2){
                               out <- lapply(tmp, as.character)[[1]]
                             }else{
                            #   last <- length(tmp)
                            #   out <- lapply(tmp, as.character)[[last]]
                               out <- lapply(tmp, as.character)
                             }
                             
                             out
                           },
                           
                           gbm =
                           {
                             #library(gbm)
                             gbmProb <- predict(modelFit, newdata, type = "response",
                                                n.trees = modelFit$tuneValue$.n.trees)
                             gbmProb[is.nan(gbmProb)] <- NA
                             
                             # need a check if all NA
                             # if so, n.trees are way too high
                             
                             
                             if(modelFit$distribution$name != "multinomial")
                             {
                               out <- ifelse(gbmProb >= .5, modelFit$obsLevels[1], modelFit$obsLevels[2])
                               ## to correspond to gbmClasses definition above
                             } else {
                               out <- colnames(gbmProb)[apply(gbmProb, 1, which.max)]
                             }
                             
                             # if there is a parameter that multiple models can be drawn, 
                             # extract these other 'lower' models
                             if(!is.null(param))
                               {
                                 tmp <- predict(modelFit, newdata, type = "response", n.trees = param$.n.trees)
                                 
                                 if(modelFit$distribution$name != "multinomial"){
                                   # if only one other parameter, need to convert to matrix
                                   if(is.vector(tmp)) tmp <- matrix(tmp, ncol = 1)
                                   tmp <- apply(tmp, 2,
                                                function(x, nm = modelFit$obsLevels) ifelse(x >= .5, nm[1], nm[2]))
                                 }else{
                                   tmp <- apply(tmp, 3,
                                                function(y, nm = modelFit$obsLevels) nm[apply(y, 1, which.max)])
                                 }
                                 
                                 # convert to list compatible splits
                                 if(length(tmp) > 1){
                                   if(!is.list(tmp)) tmp <- split(tmp, rep(1:ncol(tmp), each = nrow(tmp)))
                                 }
                                 out <- c(list(out), tmp)
                               }
                             out
                           },
                           
                           rf =
                           {
                             #library(randomForest)
                             out <-  as.character(predict(modelFit, newdata))
                             out
                           },
                           
                           svm =                           
                           {
                             #library(e1071)                             
                             out <- as.character(predict(modelFit, newdata = newdata))
                             out
                           },
                      
                           pam =
                           {
                             #library(pamr)
                             out <- as.character(
                                                 pamr.predict(modelFit,
                                                              t(newdata),
                                                              threshold = modelFit$tuneValue$.threshold))
                             #pamr.predict
                             if(!is.null(param))
                               {
                                 tmp <- vector(mode = "list", length = nrow(param) + 1)
                                 tmp[[1]] <- out
                                 for(j in seq(along = param$.threshold))
                                   {
                                     tmp[[j+1]] <- as.character(
                                                                pamr.predict(
                                                                             modelFit,
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
                               out <- predict(modelFit, newdata, s = param$.lambda, type = "class")
                               out <- as.list(as.data.frame(out, stringsAsFactors = FALSE))
                               } else {
                                 
                                 if(is.null(modelFit$lambdaOpt))
                                   stop("optimal lambda not saved; needs a single lambda value")
                                 
                                 out <- predict(modelFit, newdata, s = modelFit$lambdaOpt, type = "class")[,1]
                               }
                             out
                           },
                        
                           )
  predictedValue
}


