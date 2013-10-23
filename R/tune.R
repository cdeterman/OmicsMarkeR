
#' @title Model Tuning and Metrics
#' @description Optimizes each model based upon the parameters provided either by the internal \code{\link{denovo.grid}}
#' function or by the user.
#' @param trainVars Data used to fit the model
#' @param trainGroup Group identifiers for the training data
#' @param method A vector of strings listing models to be optimized 
#' @param k.folds Number of folds generated during cross-validation.  Default \code{"k.folds = 10"}
#' @param repeats Number of times cross-validation repeated.  Default \code{"repeats = 3"}
#' @param res Resolution of model optimization grid.  Default \code{"res = 3"}
#' @param grid Optional list of grids containing parameters to optimize for each algorithm.  
#' Default \code{"grid = NULL"} lets function create grid determined by \code{"res"}
#' @param metric Criteria for model optimization.  Available options are \code{"Accuracy"} (Predication Accuracy),
#' \code{"Kappa"} (Kappa Statistic), and \code{"AUC-ROC"} (Area Under the Curve - Receiver Operator Curve)
#' @param savePerformanceMetrics Logical argument dictating if should save the prediction data.  Default \code{savePerformanceMetrics = FALSE}
#' @param allowParallel Logical argument dictating if parallel processing is allowed via foreach package
#' @param verbose Logical argument if should output progress
#' @param theDots List of additional arguments provided in the initial classification and features selection function
#' @return Basically a list with the following elements:
#' @return \item{method}{Vector of strings listing models that were optimized}
#' @return \item{performance}{Performance generated internally to optimize model}
#' @return \item{bestTune}{List of paramaters chosen for each model}
#' @return \item{dots}{List of extra arguments initially provided}
#' @return \item{metric}{Criteria that was used for model optimization}
#' @return \item{finalModels}{The fitted models with the 'optimum' parameters}
#' @return \item{performance.metrics}{The performance metrics calculated internally for each resulting prediction}
#' @return \item{tune.metrics}{The results from each tune}
#' @return \item{perfNames}{The names of the performance metrics}
#' @return \item{comp.catch}{If the optimal PLSDA model contains only 1 component, the model must be refit with 
#' 2 components.  This catches the 1 component parameter so feature selection and further performance
#' analysis can be conducted on the 1 component.}
#' @author Charles E. Determan Jr.
# @import DiscriMiner
#' @import randomForest
#' @import caret
#' @import e1071
#' @import gbm
#' @import pamr
#' @import glmnet
#' @export

tune <- function(
  trainVars,
  trainGroup,
  method,
  k.folds = 10,
  repeats = 3,
  res = 3, 
  grid = NULL,
  metric = "Accuracy",
  savePerformanceMetrics = NULL,
  allowParallel = FALSE,
  verbose = FALSE,
  theDots = NULL)
{
  classLevels <- levels(trainGroup)
  nr <- nrow(trainVars)
  
  # for repeated cross-validation
  # creates a list of samples used for models
  repeat.index <- paste("Rep", seq(repeats), sep = "")
  for(j in 1:repeats){
    tmp <- createFolds(trainGroup, k = k.folds, list = TRUE, returnTrain = TRUE)
    names(tmp) <- paste("Fold",
                        gsub(" ", "0", format(seq(tmp))),
                        ".",
                        repeat.index[j],
                        sep = "")
    if(j == 1){
      inTrain <- tmp
      }else{
        inTrain <- c(inTrain, tmp)
      }
  }
  
  # get the remaining samples for testing group
  outTrain <- lapply(inTrain, function(inTrain, total) total[-unique(inTrain)],
                     total = seq(nr))
  
  ## combine variables and classes to make following functions simpler
  trainData <- as.data.frame(trainVars)
  trainData$.classes <- trainGroup
  
  if(is.null(grid)){
    grid <- denovo.grid(data = trainData, method = method, res = res)
  }#else{
    #verify user provided grid is okay
  #}
  
  ## get instructions to guide the loops within the modelTuner
  tune.guide <- tune.instructions(method, grid)
  
  ## run some data thru the summary function and see what we get  
  ## get phoney performance to obtain the names of the outputs
  testOutput <- vector("list", length(method))
  for(i in seq(along = method)){
    tmp <- data.frame(pred = sample(trainGroup, min(10, length(trainGroup))),
                      obs = sample(trainGroup, min(10, length(trainGroup))))
    testOutput[[i]] <- tmp
  }
  
  perfNames <- lapply(lapply(testOutput,
                             FUN = perf.calc,
                             lev = classLevels,
                             model = method), names)
  names(perfNames) <- method
  
  # tune each model  
  tmp <- modelTuner(trainData = trainData, 
                    guide = tune.guide, 
                    method = method,
                    inTrain = inTrain, 
                    outTrain = outTrain, 
                    lev = classLevels, 
                    verbose = verbose,
                    allowParallel = allowParallel,
                    theDots = theDots)

  performance <- vector("list", length(method))
  tune.results <- vector("list", length(method))
  for(i in seq(along = method)){
    performance[[i]] <- tmp[[i]]$performance
    tune.results[[i]] <- tmp[[i]]$tunes
  }
  
  tune.cm <- vector("list", length(method))
  for(i in seq(along = tune.results)){
    if(length(grep("^\\cell", colnames(tune.results[[i]]))) > 0)
    {
      tune.cm[[i]] <- tune.results[[i]][, !(names(tune.results[[i]]) %in% perfNames[[i]])]
      tune.results[[i]] <- tune.results[[i]][, -grep("^\\cell", colnames(tune.results[[i]]))]
    } else tune.cm <- NULL
    if(!is.null(tune.cm)){
      names(tune.cm) <- method
    }
  }
    
  # all possible parameter names
  paramNames <- levels(tune.guide[[1]]$model$parameter)
  
  cat("Aggregating results\n")
  flush.console()
    
  perfCols <- sapply(performance, names)
  perfCols <- lapply(perfCols, paramNames, FUN = function(x,y) x[!(x %in% y)])
  
  ## Sort the tuning parameters from least complex to most complex  
  ## lapply only works if one method being run
  #performance <- lapply(performance, "byComplexity", model = method)
  
  ## mapply only works if multiple methods being run
  #mapply("byComplexity", performance, method)
  
  #modified from caret:::byComplexity
  for(i in seq(along=method)){
    performance[[i]] <- byComplexity2(performance[[i]], method[i])
  }
  
  #if(any(is.na(performance[, metric])))
  #{
  #  warning("missing values found in aggregated results")
  #  print(performance)
  #}
  
  cat("Selecting tuning parameters\n")
  flush.console()
    
  ## select the optimal set
  #selectClass <- class(trControl$selectionFunction)[1]
  
  ## Select the "optimal" tuning parameter.
  # bestIter <- best(x = performance, metric = metric, maximize = maximize)
  #bestIter <- mapply("best", performance, metric)
  bestIter <- vector("list", length(performance))
  for(i in seq(along=performance)){
    bestIter[[i]] <- caret:::best(performance[[i]], metric, maximize=TRUE)
  }
  
  # make sure a model was chosen for each method and that it is only one option
  if(any(unlist(lapply(bestIter, is.na))) || any(unlist(lapply(bestIter, length)) != 1)) stop("final tuning parameters could not be determined")
  
  # extract the tune parameters for each model
  # added unlist to tuning parameters after modifying first lapply to choose parameter specifically
  tune.parameters <- lapply(lapply(tune.guide, "[[", 4), function(x) x["parameter"])
  bestTune <- vector("list", length(method))
  for(i in seq(along = method)){
    bestTune[[i]] <- performance[[i]][bestIter[[i]], as.character(unlist(tune.parameters[[i]])), drop = FALSE]
  }
  
  ## Save some or all of the resampling summary metrics
  # byResample <- 
  perfMetrics <- 
    if(!is.null(savePerformanceMetrics)){
      if(savePerformanceMetrics == "all"){
        # save all model metrics
        out <- tune.results
        colnames(out) <- gsub("^\\.", "", colnames(out))
        out
      }else{
        # just save the best model metrics
        out <- merge(bestTune, tune.results)
        out <- out[,!(names(out) %in% names(grid))]
        out
      }
    }else{
      out <- NULL
      out
    }
  
  ## Rename parameters to have '.' at the start of each
  #names(bestTune) <- paste(".", names(bestTune), sep = "")   
  #names(bestTune[[2]]) <- paste(".", names(bestTune[[2]]), sep="")
  newnames <- lapply(bestTune, FUN = function(x) names(x) = paste(".", names(x), sep=""))
  for(i in seq(along = newnames)){
    names(bestTune[[i]]) <- newnames[[i]]
  }
  
  ## Restore original order of performance
  for(m in seq(along = method)){
    orderList <- list()
    for(i in seq(along = tune.guide[[m]]$model$parameter))
    {
      orderList[[i]] <- performance[[m]][,as.character(tune.guide[[m]]$model$parameter[i])]
    }
    names(orderList) <- as.character(tune.guide[[m]]$model$parameter)    
    performance[[m]] <- performance[[m]][do.call("order", orderList),]
  }
  
  if(verbose)
  {
    for (i in seq(along=method)){
      cat("Fitting",
          paste(paste(gsub("^\\.", "",
                           names(bestTune[[i]])), "=",
                      bestTune[[i]]),
                collapse = ", "),
          paste("on full training set for", method[i])
          ,"\n")
    }
    flush.console()
  }
  
  # a check for plsda models - if ncomp = 1 we need to retain it to use later during feature selection as a warning
  if(any(method == "plsda")){
    catch <- which(method == "plsda")
    if(bestTune[[catch]] == 1){
      plsda.comp.catch  <- 1
    }else{
      plsda.comp.catch <- NULL
    }
  }else{
    plsda.comp.catch <- NULL
  }
  
  finalModel <- vector("list", length(method))
  for(i in seq(along = method)){
    finalModel[[i]] <- training(data = trainData,
                                method = method[i],
                                tuneValue = bestTune[[i]],
                                obsLevels = classLevels,
                                theDots = theDots)  
  }
  
  finalModel <- lapply(finalModel, function(x) x = x$fit)
  names(finalModel) <- method
  
  ## To use predict.train and automatically use the optimal lambda,
  ## we need to save it
  if(any(method %in% "glmnet")){
    m <- which(method == "glmnet")
    finalModel[[m]]$lambdaOpt <- bestTune[[m]]$.lambda
  } 
  #endTime <- proc.time()
  #times <- list(everything = endTime - startTime,
  #              final = finalTime)
  
  out <-list(
    methods = method,
    #modelType = modelType,
    performance = performance,
    #pred = tmp$predictions,
    bestTune = bestTune,
    #dots = list(...),
    dots = theDots,
    metric = metric,
    finalModels = finalModel,
    #trainingData = outData,
    #resample = perfMetrics,
    performance.metrics = perfMetrics,
    tune.metrics = tune.cm,
    perfNames = perfNames,
    #yLimits = if(is.numeric(trainGroup)) range(trainGroup) else NULL,
    comp.catch = plsda.comp.catch
    #times = times - if want time took for tuning
  )
  
  out  
}

#expand.grid(interaction.depth = seq(1, 3),
#            n.tress = floor((1:3) * 50),
#            shrinkage = c(.1, .01, .001))

