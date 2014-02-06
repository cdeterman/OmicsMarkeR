
#' @title Model Tuner
#' @description Optimizes each model based upon the parameters provided either by the internal \code{\link{denovo.grid}}
#' function or by the user.
#' @param trainData Data used to fit the model
#' @param guide Output from \code{\link{tune.instructions}}.  Facilitates the optimization by avoiding redundant model fitting.
#' @param method Vector of strins listing models to be fit
#' @param inTrain Indicies for cross-validated training models
#' @param outTrain Indicies for cross-validated testing models
#' @param lev Group levels
#' @param savePredictions Logical argument dictating if should save the prediction data.  Default \code{savePredictions = FALSE}
#' @param allowParallel Logical argument dictating if parallel processing is allowed via foreach package
#' @param verbose Logical argument dictating if should print progress
#' @param theDots List of additional arguments provided in the initial classification and features selection function
#' @return Returns list of fitted models
#' @import DiscriMiner
#' @import randomForest
#' @import plyr
#' @import caret
#' @import e1071
#' @import gbm
#' @import pamr
#' @import glmnet
#' @import foreach
#' @import caTools
# ' @export
#' @author Charles E. Determan Jr.

modelTuner <- function(trainData,
                       guide,
                       method,
                       inTrain,
                       outTrain,
                       lev,
                       savePredictions = FALSE,
                       allowParallel = FALSE,
                       verbose = FALSE,
                       theDots = NULL
)
{
  ### Parallel Programming for Windows
  #getDoParWorkers()
  #detectCores()
  #cl <- makeCluster(1)
  #registerDoSNOW(cl)
  #stopCluster(cl)
  
  `%op%` <- if(allowParallel){
    `%dopar%` 
  }else{
    `%do%`
  }  
  # set number of iterations to number of resample lists (i.e. folds * repeats)
  # iter = seq(along = inTrain)
  
  # This is a triple loop
  # first loop cycles through the methods chosen
  # second loop cycles through the cross validations
  # third loop cycles through parameters of each cross validation cycle
  #result <- foreach(iter = seq(along = inTrain), 
  #                 .combine = "c", 
  #                .verbose = FALSE, 
  #               .packages = c("methods", "caret"), 
  #              .errorhandling = "stop") %:%
  #foreach(parms = 1:nrow(guide[[i]]$loop), # how many combinations of parameters to try in these loops
  #       .combine = "c", 
  #      .verbose = FALSE, 
  #     .packages = c("methods", "caret"), 
  #    .errorhandling = "stop")  %op%
  
  #algo <- 1
  #iter <- 1
  #parms <- 1
  
  tmp.list <- foreach(algo = seq(along = method),
                      .verbose = F, 
                      .packages = c("OmicsMarkeR","foreach","caret","plyr","DiscriMiner","randomForest","e1071","gbm","pamr","glmnet","caTools"),
                      .export = c("progress", "training", "round.multiple", "predicting"),
                      .errorhandling = "stop") %:%
    foreach(iter = seq(along = inTrain), 
            .combine = "c", 
            .verbose = F, 
            .packages = c("OmicsMarkeR","foreach","caret","plyr","DiscriMiner","randomForest","e1071","gbm","pamr","glmnet","caTools"),
            .export = c("progress", "training", "round.multiple", "predicting"),
            .errorhandling = "stop") %:%
    foreach(parms = seq(nrow(guide[[algo]]$loop)), # how many combinations of parameters to try in these loops
            .combine = "c", 
            .verbose = F, 
            .packages = c("OmicsMarkeR","foreach","caret","plyr","DiscriMiner","randomForest","e1071","gbm","pamr","glmnet","caTools"),
            .export = c("progress", "training", "round.multiple", "predicting"),
            .errorhandling = "stop") %op%
{        
  ## Removes '.' from start of each parameter
  ## create 'printed' for verbose printing
  printed <- format(guide[[algo]]$loop, digits = 4)
  colnames(printed) <- gsub("^\\.", "", colnames(printed))
  
  # show progress through interations
  #caret::progress
  if(verbose) progress(printed[parms,,drop = FALSE],
                       names(inTrain), iter)
  #index <- inTrain[[iter]]
  outIndex <- outTrain[[iter]]
  
  # create models
  mod <- try(
    training(data = trainData[complete.cases(trainData[inTrain[[iter]],,drop = FALSE]),,drop = FALSE],
             method = method[algo],
             tuneValue = guide[[algo]]$loop[parms,,drop = FALSE],
             obsLevels = lev,
             theDots = theDots
    ),
    silent = TRUE)
  
  # calculate predictions if model fit successfully
  if(class(mod)[1] != "try-error")
  {          
    predicted <- try(
      predicting(method = method[algo],
                 modelFit = mod$fit,
                 orig.data = trainData,
                 indicies = inTrain[[iter]],
                 newdata = trainData[outIndex, !(names(trainData) %in% ".classes"), drop = FALSE],
                 param = guide[[algo]]$seqParam[[parms]]),
      silent = TRUE)
    
    # what should it do if predictionFunction fails???
    if(class(predicted)[1] == "try-error")
    {
      stop("prediction results failed")
      #stop(paste("prediction results failed on", method[algo], outIndex, sep = " "))
    }
    
    # what should it do if createModel fails???
  } else {
    stop("model could not be fit")
  }
  
  # If the model was built with parameters that 'submodels' can be extracted
  # this section will combine them together
  if(!is.null(guide[[algo]]$seq))
  {
    ## merge the fixed and seq parameter values together
    #caret::expandParameters
    allParam <- caret:::expandParameters(guide[[algo]]$loop[parms,,drop = FALSE], guide[[algo]]$seqParam[[parms]])
    
    ## For glmnet, we fit all the lambdas but x$fixed has an NA
    if(method[algo] == "glmnet") allParam <- allParam[complete.cases(allParam),, drop = FALSE]
    
    ## collate the predicitons across all the sub-models
    predicted <- lapply(predicted,
                        function(x, y, lv)
                        {
                          if(!is.factor(x) & is.character(x)) x <- factor(as.character(x), levels = lv)
                          data.frame(pred = x, obs = y, stringsAsFactors = FALSE)
                        },
                        y = trainData$.classes[outIndex],
                        lv = lev)
    #trainData$.classes[index]
    
    ## get the performance for this resample for each sub-model
    perf.metrics <- lapply(predicted,
                           perf.calc,
                           lev = lev,
                           model = method[algo])
    
    ## for classification, add the cell counts
    #library(plyr)
    if(length(lev) > 1)
    {
      #caret::flatTable
      cells <- lapply(predicted,
                      function(x) caret:::flatTable(x$pred, x$obs))
      for(ind in seq(along = cells)){
        perf.metrics[[ind]] <- c(perf.metrics[[ind]], cells[[ind]])
      } 
    }
    perf.metrics <- do.call("rbind", perf.metrics)          
    perf.metrics <- cbind(allParam, perf.metrics)
    
  } else {       
    
    # for models without retaining 'lower' parameters
    if(is.factor(trainData$.classes)) predicted <- factor(as.character(predicted),
                                                          levels = lev)
    tmp <-  data.frame(pred = predicted,
                       obs = trainData$.classes[outIndex],
                       stringsAsFactors = FALSE)
    
    ## Sometimes the code above does not coerce the first
    ## columnn to be named "pred" so force it
    names(tmp)[1] <- "pred"
    
    #if(ctrl$savePredictions)
    #{
    #  tmpPred <- tmp
    #  tmpPred$rowIndex <- outIndex
    #  tmpPred <- cbind(tmpPred, guide$loop[parms,,drop = FALSE])
    #  tmpPred$Resample <- names(inTrain)[iter]
    #} else tmpPred <- NULL
    
    ##################################
    
    ## Generate performance metrics
    perf.metrics <- perf.calc(tmp,
                              lev = lev,
                              model = method[algo])
    
    ## Generate the confusion matrix
    perf.metrics <- c(perf.metrics, caret:::flatTable(tmp$pred, tmp$obs))
    perf.metrics <- as.data.frame(t(perf.metrics))
    perf.metrics <- cbind(perf.metrics, guide[[algo]]$loop[parms,,drop = FALSE])
    
  }
  #perf.metrics$Resamples <- names(inTrain)[iter]
  perf.metrics$sampleIndex <- names(inTrain)[iter]
  
  # was verboseIter
  #caret::progress
  if(verbose) progress(printed[parms,,drop = FALSE],
                       names(inTrain), iter, FALSE)
  list(tunes=perf.metrics)
}
  
  ####################################
  ###### Tuning loops Complete #######
  ####################################
  if(length(method) > 1){
    names(tmp.list) <- method
    ## plyr:::rbind.fill - binds list of dataframes together
    tunes <- lapply(tmp.list, FUN = function(x) rbind.fill(x[names(x) == "tunes"]))
    #pred <- if(savePredictions)  rbind.fill(result[names(result) == "pred"]) else NULL
    
    ## remove '.' from each name
    new.names <- lapply(tunes, FUN = function(x) gsub("^\\.", "", names(x)))
    tunes <- mapply(tunes, FUN = function(x,y) {names(x) <- y; return(x)}, y = new.names, SIMPLIFY = F)
    
    for(i in length(tunes)){
      if(any(!complete.cases(tunes[[i]][,!grepl("^cell|sampleIndex", names(tunes[[i]])),drop = FALSE])))
      {
        warning(paste("There were missing values in resampled performance measures in", names(tunes[i]), sep = " "))
      } 
    }
    
    ## plyr:::ddply
    ## caret:::MeanSD
    metrics <-mapply(tunes, FUN = function(x,y) ddply(x[,!grepl("^cell|sampleIndex", colnames(x)),drop = FALSE],
                                                      y$model$parameter,
                                                      caret:::MeanSD, exclude = y$model$parameter), y = guide, SIMPLIFY = F)
    
    
    
    print("Model Tuning Complete")
    
    out <- vector("list", length(method))
    names(out) <- method
    for(i in seq(along = method)){
      out[[i]] <- list(performance = metrics[[i]], tunes = tunes[[i]])
    }
  }else{
    tmp.list <- unlist(tmp.list, recurs = FALSE)
    ## plyr:::rbind.fill - binds list of dataframes together
    tunes <- rbind.fill(tmp.list[names(tmp.list) == "tunes"])
    #pred <- if(savePredictions)  rbind.fill(result[names(result) == "pred"]) else NULL
    
    ## remove '.' from each name
    names(tunes) <- gsub("^\\.", "", names(tunes))  
    
    if(any(!complete.cases(tunes[,!grepl("^cell|sampleIndex", colnames(tunes)),drop = FALSE])))
    {
      warning("There were missing values in resampled performance measures.")
    }
    
    ## plyr:::ddply
    ## caret:::MeanSD
    metrics <- ddply(tunes[,!grepl("^cell|sampleIndex", colnames(tunes)),drop = FALSE],
                     guide[[1]]$model$parameter,
                     caret:::MeanSD, exclude = guide[[1]]$model$parameter)
    
    print(paste(method, "complete"))
    
    out <- list(performance = metrics, tunes = tunes)
  }
  
  out
}

