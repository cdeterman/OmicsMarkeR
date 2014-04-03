
#' @title Bagging Wrapper for Ensemble Features Selection
#' @description Compiles results of ensemble feature selection
#' @param X A matrix containing numeric values of each feature
#' @param Y A factor vector containing group membership of samples
#' @param method A vector listing models to be fit
#' @param bags Number of bags to be run
#' @param f Number of features desired
#' @param aggregation.metric string indicating the type of ensemble aggregation.
#' Avialable options are \code{"CLA"} (Complete Linear),
#'  \code{"EM"} (Ensemble Mean), \code{"ES"} (Ensemble Stability), and
#'  \code{"EE"} (Ensemble Exponential)
#' @param k.folds Number of folds generated during cross-validation
#' @param repeats Number of times cross-validation repeated
#' @param res Optional - Resolution of model optimization grid
#' @param tuning.grid Optional list of grids containing parameters to optimize for each algorithm.  
#' Default \code{"tuning.grid = NULL"} lets function create grid determined by \code{"res"}
#' @param optimize Logical argument determining if each model should be optimized.
#' Default \code{"optimize = TRUE"}
#' @param optimize.resample Logical argument determining if each resample should be re-optimized.
#' Default \code{"optimize.resample = FALSE"} - Only one optimization run, subsequent models use initially
#' determined parameters
#' @param metric Criteria for model optimization.  Available options are \code{"Accuracy"} (Predication Accuracy),
#' \code{"Kappa"} (Kappa Statistic), and \code{"AUC-ROC"} (Area Under the Curve - Receiver Operator Curve)
#' @param model.features Logical argument if should have number of features selected to be determined
#' by the individual model runs.  Default \code{"model.features = FALSE"}
#' @param allowParallel Logical argument dictating if parallel processing is allowed via foreach package.
#' Default \code{allowParallel = FALSE}
#' @param verbose Logical argument if should output progress
#' @param theDots Optional arguments provided for specific models or user defined parameters if 
#' \code{"optimize = FALSE"}.
#' @return \item{results}{List with the
#' following elements:}
#' @return \itemize{
#'  \item{Methods: Vector of models fit to data}
#'  \item{ensemble.results: List of length = length(method) containing aggregated features}
#'  \item{Number.bags: Number of bagging iterations}
#'  \item{Agg.metric: Aggregation method applied}
#'  \item{Number.features: Number of user-defined features}}
#' @return \item{bestTunes}{If \code{"optimize.resample = TRUE"} then returns list of 
#' best parameters for each iteration}
#' @author Charles Determan Jr
#' @import DiscriMiner
#' @import randomForest
#' @import plyr
#' @import e1071
#' @import gbm
#' @import pamr
#' @import glmnet
#' @export



bagging.wrapper <- function(X, 
                            Y, 
                            method, 
                            bags, 
                            f, 
                            aggregation.metric,
                            k.folds, 
                            repeats, 
                            res, 
                            tuning.grid,
                            optimize,
                            optimize.resample,
                            metric,
                            model.features,
                            allowParallel,
                            verbose,
                            theDots)
{
  rownames(X) <- NULL
  var.names <- colnames(X)
  nr <- nrow(X)
  nc <- ncol(X)
  # number of groups
  num.group = nlevels(Y)
  # group levels
  grp.levs <- levels(Y)  
  # leave out samples
  
  # need to retain for SVM and PAM feature selection
  trainVars.list <- vector("list", bags)
  trainGroup.list <- vector("list", bags)  
  
  if(optimize == TRUE & optimize.resample == TRUE){
    resample.tunes <- vector("list", bags)
    names(resample.tunes) <- paste("Bag", 1:bags, sep=".")
  }else{
    resample.tunes <- NULL
  }
  
  ###############
  ### Parallel Processing??
  ###############  
  for (i in 1:bags){
    # bootstrapping (i.e. random sample with replacement)
    boot=sample(nr,nr,replace=TRUE)
    
    trainVars <- X[boot,]
    trainGroup <- Y[boot]
    trainVars.list[[i]] <- trainVars
    trainGroup.list[[i]] <- trainGroup
    trainData <- as.data.frame(trainVars)
    trainData$.classes <- trainGroup
    # duplicate rownames because of bagging, must reset to 1:nrow
    rownames(trainData) <- NULL
    
    ## Run respective algorithm on bootstrapped subsamples
    if(optimize == TRUE){
      if(optimize.resample == TRUE){
        # tune the methods
        tuned.methods <- optimize.model(trainVars = trainVars,
                              trainGroup = trainGroup,
                              method = method,
                              k.folds = k.folds,
                              repeats = repeats,
                              res = res,
                              grid = tuning.grid,
                              metric = metric,
                              allowParallel = allowParallel,
                              verbose = verbose,
                              theDots = theDots)
        if(i == 1){
          finalModel <- tuned.methods$finalModel
        }else{
          finalModel <- append(finalModel, tuned.methods$finalModel)
        }
        
        # store the best tune parameters for each iteration
        names(tuned.methods$bestTune) = method
        resample.tunes[[i]] <- tuned.methods$bestTune
        
        # end of optimize.resample loop
      }else{
        if(i == 1){
          tuned.methods <- optimize.model(trainVars = trainVars,
                                    trainGroup = trainGroup,
                                    method = method,
                                    k.folds = k.folds,
                                    repeats = repeats,
                                    res = res,
                                    grid = tuning.grid,
                                    metric = metric,
                                    allowParallel = allowParallel,
                                    verbose = verbose,
                                    theDots = theDots)
          finalModel <- tuned.methods$finalModel
          names(tuned.methods$bestTune) <- method
        }else{
          # Fit remainder of resamples with initial best parameters
          #if(i == 2){
            tmp <- vector("list", length(method))
            names(tmp) <- method
          #}
    
          for(d in seq(along = method)){
            tmp[[d]] <- training(data = trainData,
                                 method = method[d],
                                 tuneValue = tuned.methods$bestTune[[d]],
                                 obsLevels = grp.levs,
                                 theDots = theDots)$fit
          }
          finalModel <- append(finalModel, tmp)
        }
      } # end of single optimization loop
      
      # end of optimizing loops
    }else{
      #theDots <- list(ncomp = 3, n.trees = 150, interaction.depth = 2, shrinkage = .01)
      names(theDots) <- paste(".", names(theDots), sep="")
      
      # sequester appropriate parameters
      args.seq <- sequester(theDots, method)
      
      # remove arguments used from theDots - also remove '.' from each
      names(theDots) <- sub(".", "", names(theDots))
      moreDots <- theDots[!names(theDots) %in% args.seq$pnames]
      if(length(moreDots) == 0){
        moreDots <- NULL
      }
      
      #moreDots <- theDots[-names(args.seq)]
      
      finalModel <- vector("list", length(method))
      for(q in seq(along = method)){
        finalModel[[q]] <- training(data = trainData,
                                    method = method[q],
                                    tuneValue = args.seq$parameters[[q]],
                                    obsLevels = grp.levs,
                                    theDots = moreDots)  
      }
      # end of non-optimized model fitting
    }
    
  # end of bagging loop
  }
  
  # sort models together (e.g. first 5 are "plsda", next 5 "gbm", etc.)    
  method.names <- unlist(lapply(method, FUN = function(x) paste(c(rep(x, bags)), seq(bags), sep = ".")))
  #orig.method.names <- unlist(lapply(method, FUN = function(x) c(rep(x, bags))))
  
  names(finalModel) <- paste(method, rep(seq(bags), each = length(method)), sep = ".")
  finalModel <- finalModel[match(method.names, names(finalModel))]    
  #names(finalModel) <- orig.method.names
  
  # Create empty list for features identified by each chosen algorithm
  features <- vector("list", length(method))
  names(features) <- tolower(method)
  
  
  for(j in seq(along = method)){
    ### Extract important features
    # pam requires a special mydata argument
    mydata <- vector("list", bags)
    if(method[j] == "pam"){
      for(t in 1:bags){
        mydata[[t]] <- list(x=t(trainVars.list[[t]]), y=factor(trainGroup.list[[t]]), geneid = as.character(colnames(trainVars.list[[t]])))
      }
    }else{
      # svm requires training data for RFE
      for(t in 1:bags){
        mydata[[t]] <- trainVars.list[[t]]
      }
    }
    
    if(j == 1){
      start <- 1
      end <- bags
    }
    
    if(method[j] == "svm" | method[j] == "pam" | method[j] == "glmnet"){
      bt <- vector("list", bags)
      for(l in seq(bags)){
        if(optimize == TRUE){
          if(optimize.resample == FALSE){
            bt[[l]] <- tuned.methods$bestTune[[j]]
          }else{
            bt[[l]] <- tuned.methods$bestTune[[l]]
          }
        }
      }
    }else{
      bt <- vector("list", bags)
    }
    
    if(method[j] == "plsda"){
      cc <- vector("list", bags)
      for(c in seq(bags)){
        if(optimize == TRUE){
          if(optimize.resample == FALSE){
            cc[[c]] <- tuned.methods$bestTune[[j]]
          }else{
            cc[[c]] <- tuned.methods$bestTune[[c]]
          }
        }
      }
    }
    
    finalModel.bag <- finalModel[start:end]
    tmp <- vector("list", bags)
    for(s in seq(bags)){
      tmp[[s]] <- extract.features(
        x = finalModel.bag[s],
        dat = mydata[[s]],
        grp = trainGroup.list[[s]],
        # add in gbm best tune trees???
        bestTune = bt[[s]],
        model.features = FALSE, 
        method = method[j], 
        f = NULL, 
        #similarity.metric = similarity.metric,
        comp.catch = cc)
    }
    
    if(method[j] == "glmnet"){
      features[[j]] <- data.frame(do.call("cbind", unlist(unlist(tmp, recursive = F), recursive = F)))
    }else{
      features[[j]] <- do.call("cbind", unlist(tmp, recursive = F))
      if(class(features[[j]]) != "data.frame"){
        features[[j]] <- data.frame(features[[j]])
      }
    }
    rownames(features[[j]]) <- colnames(X)
    
    start <- start + bags
    end <- end + bags
  }
  
  ### Ensemble Aggregation
  #convert to numeric & set rownames
  features.num <- lapply(features, FUN = function(z) sapply(z, FUN = function(x) as.numeric(as.character(x))))
  features.num <- lapply(features.num, function(x) {
    rownames(x) <- var.names
    return(x)
  })
  
  # Generate summary lists of each algorithm
  agg <- lapply(features.num, FUN = function(x) aggregation(efs = x, metric = aggregation.metric, f = f))
  
  # Basic ensemble model parameters
  ensemble.results <- list(Methods = method,
                           ensemble.results = agg,
                           Number.Bags = bags,
                           Agg.metric = aggregation.metric,
                           Number.features = f)
  
  
  out <- list(results = ensemble.results,
              bestTunes = resample.tunes)
  out
}

