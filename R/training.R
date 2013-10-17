training <-
  function(data, method, tuneValue, obsLevels, theDots = NULL)
    {
    
    #function(variables, classes, method, tuneValue, obsLevels, pp = NULL, last = FALSE, custom = NULL, classProbs, ...)
    data <- as.data.frame(data)

  # make sure '.classes' is set as factor
    
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

  ## the outcome data better be a factor!!!
  #type <- if(is.factor(data$.outcome)) "Classification"
  
  ## We refactor the class labels. Some methods bark/crash when there are
  ## factor levels that do not have values represented in the data (nnet produces
  ## a warning and randomForest throws an error). 
  #if(type == "Classification") data$.outcome <- factor(as.character(data$.outcome), levels = obsLevels)
  data$.classes <- factor(as.character(data$.classes), levels = obsLevels)

  xNames <- names(data)[!(names(data) %in% ".classes")]
  
  ## Later, when we assign predictions, we will convert predictions to 
  ## character and then create factors from them with levels originally
  ## found in the object obsLevels.

  trainX <- data[,!(names(data) %in% ".classes"), drop = FALSE]
  trainY <- data[,".classes"]
  
  if(tolower(method) == "gbm" & length(obsLevels) == 2)  numClasses <- ifelse(data$.classes == obsLevels[1], 1, 0)
  
  #numClasses <- ifelse(trainData$.classes == obsLevels[1], 1, 0)
  
  modelFit <- switch(tolower(method),
                     
                     ### DiscriMiner:::plsDA fits model and predictions simultaneously
                     # it is skipped during training but used during final model fitting
                     plsda =
                       {
                         #require(DiscriMiner)
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
                       require(gbm)
                       ## train will figure out whether we are doing classification or reggression
                       ## from the class of the outcome and automatically specify the value of
                       ## 'distribution' in the control file. If the user wants to over-ride this,
                       ## this next bit will allow this.
                       
                       
                       # need to make sure only extract arguments that pertain to gbm
                       gbm.args <- c("w", "var.monotone", "n.minobsinnode", 
                                     "bag.fraction", "var.names", "response.name", "group") 
                       theDots <- theDots[names(theDots) %in% gbm.args]
                       
                       #theDots <- list(...)
                       #if(any(names(theDots) == "distribution"))
                         #{
                         #  modDist <- theDots$distribution
                          # theDots$distribution <- NULL
                         #} else #{
                        #   if(type == "Regression")
                         #    {
                          #     modDist <- "gaussian"
                           #  } else modDist <- if(length(obsLevels) == 2)  "bernoulli" else "multinomial"
                         #}
                       gbmdist <- if(length(unique(trainY)) == 2){
                         "bernoulli"}else{
                           "multinomial"
                         }         
                       
                       # check gbm setup file to see if this is necessary
                       modY <- if(gbmdist != "multinomial") numClasses else trainY
                       
                       if(gbmdist != "multinomial"){
                         modY <- numClasses
                         #gbm.dat <- data.frame(trainX, modY)
                       }else{
                         modY <- trainY
                         #gbm.dat <- data.frame(trainX, modY)
                       }
                       
                       modArgs <- list(x = trainX,
                                       y = modY,
                                       #data = gbm.dat,
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
                       require(randomForest)
                       
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
                       
                      # randomForest(trainX, 
                      #              trainY, 
                      #              mtry = tuneValue$.mtry, 
                      #              ntree=round.multiple(sqrt(ncol(trainX))),
                      #              ...)
                     },                  
                    
                     # SVM will have it's tuining internal
                     # Therefore it should be removed from this loop???
                     # Is this feasible or should it just be initially optimized?
                     svm =   
                       { 
                         require(e1071)
                         out <- svm(trainX, 
                                    trainY,
                                    cost = 10,
                                    #cost = tuneValue$.C, 
                                    cachesize=500,
                                    scale=F, 
                                    type="C-classification", 
                                    kernel="linear")                         
                         out
                       },
  
                     pam = 
                     {
                       require(pamr)    
                       
                       pamr.args <- c("n.threshold", "threshold.scale", "scale.sd", "se.scale")
                       theDots <- theDots[names(theDots) %in% pamr.args]
                       
                       modArgs <- list(data = list(x = t(trainX), y = trainY, geneid = as.character(colnames(trainX))),
                                       threshold = tuneValue$.threshold
                                       )
                       
                       if(length(theDots) > 0) modArgs <- c(modArgs, theDots)
                       
                       do.call("pamr.train", modArgs)
                       
                      # out <- pamr.train(list(x = t(trainX), y = trainY),
                      #                   threshold = tuneValue$.threshold, 
                      #                   ...)
                      # out
                     },         
                      
                     # need to address glmnet still
                     glmnet =
                     {
                       require(glmnet)
                       numLev <- if(is.character(trainY) | is.factor(trainY)) length(levels(trainY)) else NA

                       #theDots <- list(...)
                       
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
                       }
                       
                       modelArgs <- c(
                                      list(
                                           x = as.matrix(trainX),
                                           y = trainY,
                                           alpha = tuneValue$.alpha),
                                      theDots)
                       
                       out <- do.call("glmnet", modelArgs) 
                       out 
                     }
                     )
  
  
  ##save a few items so we have a self contained set of information in the model. this will
  ## also carry over to the finalModel if returnData = TRUE in train call

  
  ###### May not need to deal with S4 classes ############3
  
  ## for models using S4 classes, you can't easily append data, so 
  ## exclude these and we'll use other methods to get this information
  if(!(tolower(method) %in% tolower(c("svmLinear"))))
    {
      modelFit$xNames <- xNames
    #  modelFit$problemType <- type
      modelFit$tuneValue <- tuneValue
      modelFit$obsLevels <- obsLevels
    }
  #if(!is.null(modelFit) && 
  #   any(names(modelFit) == "call" & 
  #   !(method %in% c("rpart", "rpart2", "earth", "fda")))) 
  #     modelFit$call <- scrubCall(modelFit$call)
  #require(methods)
  #if(length(slotNames(modelFit)) > 0 && any(slotNames(modelFit) == "call")) modelFit@call <- scrubCall(modelFit@call)

  list(fit = modelFit)
}
