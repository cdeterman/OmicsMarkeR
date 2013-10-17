predicting <- function(method, modelFit, orig.data, indicies, newdata, param = NULL)
{
  if(any(colnames(newdata) == ".classes")) newdata$.classes <- NULL

  coerceChar <- function(x){
    as.data.frame(lapply(x, as.character), stringsAsFactors = FALSE)
  }  
  
  #if(!is.null(preProc)) newdata <- predict(preProc, newdata)
  
  predictedValue <- switch(tolower(method),
                           ### plsda removed because DiscriMiner fits model and prediction simultaneously
                           
                           plsda =
                           {
                             # require(DiscriMiner)
                             # retain.models omitted because when this is used, the final model is only using the best component
                             # may switch to retain all models but this omits more processing that is likely superfluous
                             
                             # check for number of components provided.  This is important following selection of the best model
                             if(param$.ncomp == 1){
                               warning("PLSDA model contained only 1 component. PLSDA requires at least 2 components.\nModel fit with 2 components")
                               param$.ncomp = 2
                             }                             
                             
                             tmp <- plsDA(orig.data[,-which(names(orig.data) %in% c(".classes"))], 
                                          orig.data[,c(".classes")],
                                          autosel=F,
                                          #learn = orig.data[,-which(names(orig.data) %in% c(".classes"))][inTrain,, drop = F],
                                          learn = indicies,
                                          test = seq(nrow(orig.data))[-unique(indicies)],
                                          validation = "learntest",
                                          #comps = 2,
                                          comps = param$.ncomp,
                                          cv ="none",
                                          retain.model = TRUE)$classification
                             
                             if(param$.ncomp < 2){
                               out <- lapply(tmp, as.character)[[1]]
                             }else{
                               last <- length(tmp)
                               out <- lapply(tmp, as.character)[[last]]
                             }
                             out
                           },
                           
                           gbm =
                           {
                             require(gbm)
                             gbmProb <- predict(modelFit, newdata, type = "response",
                                                n.trees = modelFit$tuneValue$.n.trees)
                             gbmProb[is.nan(gbmProb)] <- NA
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
                                 if(!is.list(tmp)) tmp <- split(tmp, rep(1:ncol(tmp), each = nrow(tmp)))
                                 out <- c(list(out), tmp)
                               }
                             out
                           },
                           
                           rf =
                           {
                             require(randomForest)
                             out <-  as.character(predict(modelFit, newdata))
                             out
                             #print(out)
                           },
                           
                           svm =                           
                           {
                             require(e1071)
                             out <- as.character(predict(modelFit, newdata))
                             out
                           },
                      
                           pam =
                           {
                             require(pamr)
                             out <- as.character(
                                                 pamr.predict(modelFit,
                                                              t(newdata),
                                                              threshold = modelFit$tuneValue$.threshold))
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
                             require(glmnet)
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


