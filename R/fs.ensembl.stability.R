
#' @title Ensemble Classification & Feature Selection
#' @description Applies ensembles of models to high-dimensional data to both classify and determine important
#' features for classification.  The function bootstraps a user-specified number of times to facilitate
#' stability metrics of features selected thereby providing an important metric for biomarker investigations,
#' namely whether the important variables can be identified if the models are refit on 'different' data.
#' @param X A matrix containing numeric values of each feature
#' @param Y A factor vector containing group membership of samples
#' @param method A vector listing models to be fit.
#' Available options are \code{"plsda"} (Partial Least Squares Discriminant Analysis),
#'  \code{"rf"} (Random Forest), \code{"gbm"} (Gradient Boosting Machine),
#'  \code{"svm"} (Support Vector Machines), \code{"glmnet"} (Elastic-net Generalized Linear Model),
#'  and \code{"pam"} (Prediction Analysis of Microarrays)
#' @param k Number of bootstrapped interations
#' @param p Percent of data to by 'trained'
#' @param f Number of features desired.  Default is top 10% 
#' \code{"f = ceiling(ncol(variables)/10)"}.
#' If rank correlation is desired, set \code{"f = NULL"}
#' @param bags Number of iterations for ensemble bagging.  Default \code{"bags = 40"}
#' @param aggregation.metric String indicating which aggregation metric for features selected during bagging.
#' #' Avialable options are \code{"CLA"} (Complete Linear),
#'  \code{"EM"} (Ensemble Mean), \code{"ES"} (Ensemble Stability), and
#'  \code{"EE"} (Ensemble Exponential)
#' @param stability.metric string indicating the type of stability metric.
#' Avialable options are \code{"jaccard"} (Jaccard Index/Tanimoto Distance),
#'  \code{"sorensen"} (Dice-Sorensen's Index), \code{"ochiai"} (Ochiai's Index),
#'  \code{"pof"} (Percent of Overlapping Features), \code{"kuncheva"} (Kuncheva's Stability Measures),
#'  \code{"spearman"} (Spearman Rank Correlation), and \code{"canberra"} (Canberra Distance)
#' @param optimize Logical argument determining if each model should be optimized.
#' Default \code{"optimize = TRUE"}
#' @param optimize.resample Logical argument determining if each resample should be re-optimized.
#' Default \code{"optimize.resample = FALSE"} - Only one optimization run, subsequent models use initially
#' determined parameters
#' @param tuning.grid Optional list of grids containing parameters to optimize for each algorithm.  
#' Default \code{"tuning.grid = NULL"} lets function create grid determined by \code{"res"}
#' @param k.folds Number of folds generated during cross-validation.  Default \code{"k.folds = 10"}
#' @param repeats Number of times cross-validation repeated.  Default \code{"repeats = 3"}
#' @param resolution Optional - Resolution of model optimization grid.  Default \code{"res = 3"}
#' @param metric Criteria for model optimization.  Available options are \code{"Accuracy"} (Predication Accuracy),
#' \code{"Kappa"} (Kappa Statistic), and \code{"AUC-ROC"} (Area Under the Curve - Receiver Operator Curve)
#' @param model.features Logical argument if should have number of features selected to be determined
#' by the individual model runs.  Default \code{"model.features = FALSE"}
#' @param allowParallel Logical argument dictating if parallel processing is allowed via foreach package.
#' Default \code{allowParallel = FALSE}
#' @param verbose Logical argument if should output progress
#' @param ... Extra arguments that the user would like to apply to the models
#' 
#' @return \item{Methods}{Vector of models fit to data}
#' @return \item{performance}{Performance metrics of each model and bootstrap iteration}
#' @return \item{RPT}{Robustness-Performance Trade-Off as defined in Saeys 2008}
#' @return \item{features}{List concerning features determined via each algorithms feature selection criteria.}
#' @return \itemize{
#'  \item{metric: Stability metric applied}
#'  \item{features: Matrix of selected features}
#'  \item{stability: Matrix of pairwise comparions and average stability}
#'  }
#' @return \item{stability.models}{Function perturbation metric - i.e. how similar are the features selected
#' by each model.}
#' @return \item{all.tunes}{If \code{"optimize.resample = TRUE"} then returns list of 
#' optimized parameters for each bagging and bootstrap interation.}
#' @return \item{final.best.tunes}{If \code{"optimize.resample = TRUE"} then returns list of
#' optimized parameters for each bootstrap of the bagged models refit to aggregated selected features.}
#' @return \item{specs}{List with the
#' following elements:}
#' @return \itemize{
#'  \item{total.samples: Number of samples in original dataset}
#'  \item{number.features: Number of features in orginal dataset}
#'  \item{number.groups: Number of groups}
#'  \item{group.levels: The specific levels of the groups}
#'  \item{number.observations.group: Number of observations in each group}}
#' @author Charles Determan Jr
#' @references Saeys Y., Abeel T., et. al. (2008) \emph{Machine Learning and Knowledge Discovery in Databases}. 
#' 313-325. http://link.springer.com/chapter/10.1007/978-3-540-87481-2_21
#' @import DiscriMiner
#' @import randomForest
#' @import plyr
#' @import caret
#' @import e1071
#' @import gbm
#' @import pamr
#' @import glmnet
#' @export

fs.ensembl.stability <- 
  function(X,
           Y,
           method,
           k = 10,
           p = 0.9,
           f = ceiling(ncol(X)/10),
           bags=40,
           aggregation.metric = "CLA",
           stability.metric = "jaccard",
           optimize = TRUE,
           optimize.resample = FALSE,
           tuning.grid = NULL,
           k.folds = if(optimize) 10 else NULL,
           repeats = if(optimize) 3 else NULL,
           resolution = if(optimize) 3 else NULL,
           metric = "Accuracy",
           model.features = FALSE,
           allowParallel = FALSE,
           verbose = FALSE,
           ...
           )
    {
    verify_data <- verify(x = X, y = Y, method = method, f = f, stability.metric = stability.metric, model.features = model.features, na.rm = FALSE)
    #verify_data <- my_verify(variables, groups, na.rm = FALSE)
    
    if(model.features == TRUE){
      stop("Error... Model derived features cannot be used for ensemble because all features are ranked and then subset.\nSet model.features = FALSE")
    }
    
    X <- verify_data$X
    Y <- verify_data$Y
    method <- verify_data$method
    f <- verify_data$f
    
    raw.data <- as.data.frame(X)
    raw.data$.classes <- Y
    
    method <- tolower(method)
    nr<-nrow(X)
    nc<-ncol(X)
    # number of groups
    num.group <- nlevels(Y)
    # what the groups are
    grp.levs <- levels(Y)  
    # how many obs in each group
    num.obs.group <- as.vector(table(Y))
    theDots <- list(...)
    
    # Create empty list for features identified by each chosen algorithm
    features <- vector("list", k)
    names(features) <- paste("Resample", seq(k), sep = ".")
    
    if(optimize == TRUE & optimize.resample == TRUE){
      bagged.tunes <- vector("list", k)
      names(bagged.tunes) <- paste("Resample", seq(k), sep=".")
    }else{
      bagged.tunes <- NULL
    }    
    
    inTrain <- rlply(k, createDataPartition(Y, p = p, list = FALSE))
    #inTrain <- rlply(k, sample(nr, round(p*nr)))
    outTrain <- lapply(inTrain, function(inTrain, total) total[-unique(inTrain)],
                       total = seq(nr))    
    
    ### Stability Loop
    for (i in 1:k){      
      trainX <- X[inTrain[[i]],, drop=F]
      trainY <- Y[inTrain[[i]], drop=F]
      trainData <- as.data.frame(trainX)
      trainData$.classes <- trainY
      
      ## Bagging loop
      results.bagging <- bagging.wrapper(X = trainX, 
                                         Y = trainY, 
                                         method = method,
                                         bags = bags, 
                                         f = f, 
                                         aggregation.metric = aggregation.metric,
                                         k.folds = k.folds,
                                         repeats = repeats,
                                         res = resolution,
                                         tuning.grid = tuning.grid,
                                         optimize = optimize,
                                         optimize.resample = optimize.resample,
                                         metric = metric,
                                         model.features = model.features,
                                         verbose = verbose,
                                         allowParallel = allowParallel,
                                         theDots = theDots)
      
      # store the best tune parameters for each iteration
      if(!is.null(bagged.tunes)){
        bagged.tunes[[i]] <- results.bagging$bestTunes
      }
      
      # Store the features selected for stability analysis
      features[[i]] <- as.list(results.bagging$results$ensemble.results)
      
      ### Re-fitting models to reduced features
      # subset only those features which were selected
      trainData.new <- lapply(features[[i]], FUN = function(x) trainData[,colnames(trainData) %in% c(as.vector(names(x)), ".classes")])
      
      if(optimize == TRUE){
        if(optimize.resample == TRUE){
          tunedModel.new <- vector("list", length(method))
          for(m in seq(along = method)){
            tunedModel.new[[m]] <- optimize(trainVars = trainData.new[[m]][,!colnames(trainData.new[[m]]) %in% c(".classes")],
                                        trainGroup = trainData.new[[m]]$.classes,
                                        method = method[m],
                                        k.folds = k.folds,
                                        repeats = repeats,
                                        res = resolution,
                                        grid = tuning.grid,
                                        metric = metric,
                                        savePerformanceMetrics = NULL,
                                        allowParallel = allowParallel,
                                        verbose = verbose,
                                        theDots = theDots)            
          }
          
          if(i == 1){
            finalModel.new <- sapply(tunedModel.new, FUN = function(x) x$finalModel)
            new.best.tunes <- sapply(tunedModel.new, FUN = function(x) x$bestTune)
            names(new.best.tunes) <- method
          }else{
            tmp.model <- sapply(tunedModel.new, FUN = function(x) x$finalModels)
            tmp.tunes <- sapply(tunedModel.new, FUN = function(x) x$bestTune)
            names(tmp.tunes) <- method
            finalModel.new <- c(finalModel.new, tmp.model)
            new.best.tunes <- append(new.best.tunes, tmp.tunes)
          }  
          # end of full optimize loop
        }else{
          if(i == 1){
            tunedModel.new <- vector("list", length(method))
            for(m in seq(along = method)){
              tunedModel.new[[m]] <- optimize(trainVars = trainData.new[[m]][,!colnames(trainData.new[[m]]) %in% c(".classes")],
                                          trainGroup = trainData.new[[m]]$.classes,
                                          method = method[m],
                                          k.folds = k.folds,
                                          repeats = repeats,
                                          res = resolution,
                                          grid = tuning.grid,
                                          metric = metric,
                                          savePerformanceMetrics = NULL,
                                          allowParallel = allowParallel,
                                          verbose = verbose,
                                          theDots = theDots)  
            
            }
          }else{
            tmp <- vector("list", length(method))
            names(tmp) <- method
            for(d in seq(along = method)){
              tmp[[d]] <- training(data = trainData.new[[d]],
                                   method = method[d],
                                   tuneValue = data.frame(t(unlist(tunedModel.new[[d]]$bestTune))),
                                   obsLevels = grp.levs,
                                   theDots = theDots)$fit
            }
          }
          
          if(i == 1){
            finalModel.new <- sapply(tunedModel.new, FUN = function(x) x$finalModels)
            new.best.tunes <- sapply(tunedModel.new, FUN = function(x) x$bestTune)
          }else{
            tmp.model <- lapply(tmp, FUN = function(x) x)
            finalModel.new <- c(finalModel.new, tmp.model)
          }  
        } # end of single optimize loop
        # end of optimize loop
      }else{
        names(theDots) <- paste(".", names(theDots), sep="")
        
        # sequester appropriate parameters to fit models
        args.seq <- sequester(theDots, method)
        
        # remove arguments used from theDots - also remove '.' from each
        names(theDots) <- sub(".", "", names(theDots))
        moreDots <- theDots[!names(theDots) %in% args.seq$pnames]
        if(length(moreDots) == 0){
          moreDots <- NULL
        }
        
        tunedModel.new <- vector("list", length(method))
        for(q in seq(along = method)){
          tunedModel.new[[q]] <- training(data = trainData.new[[q]],
                                          method = method[q],
                                          tuneValue = args.seq$parameters[[q]],
                                          obsLevels = grp.levs,
                                          theDots = moreDots)$fit  
        }
        
        if(i == 1){
          finalModel.new <- sapply(tunedModel.new, FUN = function(x) x)
        }else{
          tmp.model <- lapply(tunedModel.new, FUN = function(x) x)
          finalModel.new <- c(finalModel.new, tmp.model)
        }
      } # end of non-optimized loop
    } # end stability loop
    
    final.features <- features
    #ref.features <- lapply(features, FUN = function(x) data.frame(x))
    
    ### Performance Metrics of Reduced Models
    final.metrics <- prediction.metrics(finalModel = finalModel.new,
                                        method = method,
                                        raw.data = raw.data,
                                        inTrain = inTrain,
                                        outTrain = outTrain,
                                        bestTune = new.best.tunes,
                                        features = final.features,
                                        grp.levs = grp.levs)
    
    
    ### Extract Performance Metrics
    if(optimize == TRUE){
      if(optimize.resample == TRUE){        
        colnames(final.metrics) <- gsub("^\\.", "", colnames(final.metrics))  
        performance <- vector("list", length(method))
        names(performance) <- method
        final.metrics <- final.metrics[,!grepl("^cell", colnames(final.metrics)),drop = FALSE]
        
        # sort the list elements so applied in the proper order
        method.names <- unlist(lapply(method, FUN = function(x) c(rep(x, length(new.best.tunes)/length(method)))))
        new.best.tunes <- new.best.tunes[order(names(new.best.tunes), levels = method.names)]
        
        # dataframe retains colnames after split
        rownames(final.metrics) <- letters[seq(nrow(final.metrics))]
        test.final.metrics <- split(as.data.frame(final.metrics), rownames(final.metrics))
        for(i in seq(along = test.final.metrics)){
          rownames(test.final.metrics[[i]]) <- 1
        }
        rownames(final.metrics) <- method.names
        all.model.perfs <- mapply(new.best.tunes, FUN = function(x, y) list(x,y), y = test.final.metrics, SIMPLIFY = FALSE)
        
        # name resamples within all.model.perfs
        samps <- rep(rep(paste("Resample", seq(k), sep = "."), each = 2), length(method))
        for(i in seq(length(all.model.perfs))){
          if(i == 1){
            start = i
            end = i+1
          }
          names(all.model.perfs[[i]]) <- samps[start:end]
          start = start + 2
          end = end + 2
        }
        
        for(h in seq(along = method)){
          x <- final.metrics[,!grepl("^cell", colnames(final.metrics)),drop = FALSE]
          tmp <- subset(x, rownames(x) == method[h])
          performance[[h]] <- c(colMeans(tmp, na.rm = TRUE), apply(x, 2, sd, na.rm = TRUE))
          names(performance[[h]])[-(1:ncol(tmp))] <- paste(names(performance[[h]])[-(1:ncol(tmp))], "SD", sep = " ")
          performance[[h]] <- t(data.frame(performance[[h]]))
          rownames(performance[[h]]) <- 1
        }
      }else{
        colnames(final.metrics) <- gsub("^\\.", "", colnames(final.metrics))  
        performance <- vector("list", length(method))
        names(performance) <- method
        for(h in seq(along = method)){
          x <- final.metrics[,!grepl("^cell", colnames(final.metrics)),drop = FALSE]
          #x <- x[, !colnames(x) %in% rownames(final.metrics), drop = FALSE]
          tmp <- subset(x, rownames(x) == method[h])
          performance[[h]] <- c(colMeans(tmp, na.rm = TRUE), apply(x, 2, sd, na.rm = TRUE))
          names(performance[[h]])[-(1:ncol(tmp))] <- paste(names(performance[[h]])[-(1:ncol(tmp))], "SD", sep = " ")
          #performance[[h]] <- do.call(cbind, c(as.vector(tuned.methods$bestTune[[h]]), performance[[h]]))
          performance[[h]] <- do.call(cbind, c(as.vector(new.best.tunes[[h]]), performance[[h]]))
          colnames(performance[[h]]) <- gsub("^\\.", "", colnames(performance[[h]]))
          rownames(performance[[h]]) <- 1
        }
      }
    }else{
      colnames(final.metrics) <- gsub("^\\.", "", colnames(final.metrics))  
      performance <- vector("list", length(method))
      names(performance) <- method
      for(h in seq(along = method)){
        x <- final.metrics[,!grepl("^cell", colnames(final.metrics)),drop = FALSE]
        #x <- x[, !colnames(x) %in% rownames(final.metrics), drop = FALSE]
        tmp <- subset(x, rownames(x) == method[h])
        performance[[h]] <- c(colMeans(x, na.rm = TRUE), apply(x, 2, sd, na.rm = TRUE))
        names(performance[[h]])[-(1:ncol(tmp))] <- paste(names(performance[[h]])[-(1:ncol(tmp))], "SD", sep = " ")
        performance[[h]] <- do.call(cbind, c(as.vector(args.seq$parameters[[h]]), performance[[h]]))
        colnames(performance[[h]]) <- gsub("^\\.", "", colnames(performance[[h]]))
        rownames(performance[[h]]) <- 1
      }
    }    
    
    # need to split features into length(method) dataframes for pairwise.stability
    results.stability <- vector("list", length(method))
    names(results.stability) <- method
    for(c in seq(along = method)){
      if(stability.metric %in% c("jaccard", "kuncheva", "pof", "ochiai", "sorenson")){
        results.stability[[c]] <- as.data.frame(sapply(features, FUN = function(x) names(x[[c]])))  
      }else{
        tmp <- lapply(features, FUN = function(x) x[[c]])
        myMat <- matrix(NA, nrow = nc, ncol = k,
                        dimnames = list(colnames(X)))
        for (z in seq_along(tmp)) {
          myMat[names(tmp[[z]]), z] <- tmp[[z]]
        }
        colnames(myMat) <- paste("Resample", seq(k), sep = ".")
        results.stability[[c]] <- myMat
      }
    }
    
    # Calculate All Pairwise Stability Measurements for each algorithm's set of features
    stability <- lapply(results.stability, pairwise.stability, stability.metric = stability.metric, nc= nc)
    
    # stability across algorithms (i.e. 'function perturbation' ensemble analysis)
    if(length(method) > 1){
      stability.models <- pairwise.model.stability(features = results.stability,
                                                   stability.metric = stability.metric,
                                                   m = length(method),
                                                   k = k,
                                                   nc = nc)
    }else{
      stability.models <- NULL
    }
    
    # harmonic mean of stability and performance
    rpt.stab <- lapply(stability, FUN = function(x) x$overall)
    rpt.perf <- lapply(performance, FUN = function(x) as.data.frame(x)$Accuracy)
    rpt <- mapply(rpt.stab, FUN = function(x,y) RPT(stability = x, performance = y), y = rpt.perf)
    
    
    # add stability metric to each respective algorithm
    for(i in 1:length(method)){
      results.stability[[i]] <- list(metric = stability.metric, features = results.stability[[i]], stability = stability[[i]])
    }
    
    specs = list(total.samples=nr, 
                 number.features=nc, 
                 number.groups=num.group, 
                 group.levels=grp.levs, 
                 number.observations.group=num.obs.group)
    
    ## add remainder of data
    overall <- list(methods = method,                      # algorithms run
                    performance = performance,
                    RPT = rpt,
                    features = results.stability,          # list of each algorithms results
                    stability.models = stability.models,   # stability amongst algorithms
                    all.tunes = bagged.tunes,              # if optimize.resample returns the best tunes for each iteration
                    final.best.tunes = if(optimize.resample) all.model.perfs else NULL,   # if optimize.resample, provide parameter with performance statistics
                    specs = specs                          # general specs of the input data
    )
    return(overall)
}
