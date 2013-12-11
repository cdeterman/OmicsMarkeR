
#' @title Classification & Feature Selection
#' @description Applies models to high-dimensional data to both classify and determine important
#' features for classification.  The function bootstraps a user-specified number of times to facilitate
#' stability metrics of features selected thereby providing an important metric for biomarker investigations,
#' namely whether the important variables can be identified if the models are refit on 'different' data.
#' @param X A scaled matrix or dataframe containing numeric values of each feature
#' @param Y A factor vector containing group membership of samples
#' @param method A vector listing models to be fit.
#' Available options are \code{"plsda"} (Partial Least Squares Discriminant Analysis),
#'  \code{"rf"} (Random Forest), \code{"gbm"} (Gradient Boosting Machine),
#'  \code{"svm"} (Support Vector Machines), \code{"glmnet"} (Elastic-net Generalized Linear Model),
#'  and \code{"pam"} (Prediction Analysis of Microarrays)
#' @param k Number of bootstrapped interations
#' @param p Percent of data to by 'trained'
#' @param f Number of features desired.  Default is top 10% 
#' \code{"f = ceiling(ncol(X)/10)"}.
#' If rank correlation is desired, set \code{"f = NULL"}
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
#' @param resolution Resolution of model optimization grid.  Default \code{"resolution = 3"}
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
#' @return \item{original.best.tunes}{If \code{"optimize.resample = TRUE"} then returns list of 
#' optimized parameters for each bootstrap.}
#' @return \item{final.best.tunes}{If \code{"optimize.resample = TRUE"} then returns list of
#' optimized parameters for each bootstrap of models refit to selected features.}
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
# ' @import tcltk
#' @export

fs.stability <- 
  function(X,
           Y,
           method,
           k = 10,
           p = 0.9,
           f = ceiling(ncol(X)/10),
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
    #### Filter Methods???
    ## Bhattacharyya Distance
    #??bhattacharyya
    #require(fpc)
    #bhattacharyya.dist
    ## Relative Entropy
    
    verify_data <- verify(x = X, y = Y, method = method, f = f, stability.metric = stability.metric, model.features = model.features, na.rm = FALSE)
    
    X <- verify_data$X
    Y <- verify_data$Y
    method <- verify_data$method
    f <- verify_data$f
    
    raw.data <- as.data.frame(X)
    raw.data$.classes <- Y
    
    nr <- nrow(X)
    nc <- ncol(X)
    # number of groups
    num.group <- nlevels(Y)
    # what the groups are
    grp.levs <- levels(Y)  
    # how many obs in each group
    num.obs.group <- as.vector(table(Y))
    theDots <- list(...)
    
    # Create empty list for features identified by each chosen algorithm
    final.features <- vector("list", k)
    names(final.features) <- paste("Resample", seq(k), sep = ".")
    
    # need to retain for SVM and PAM feature selection
    trainVars.list <- vector("list", k)
    trainGroup.list <- vector("list", k)  
    
    if(optimize == TRUE & optimize.resample == TRUE){
      resample.tunes <- vector("list", k)
      names(resample.tunes) <- paste("Resample", 1:k, sep=".")
    }else{
      resample.tunes <- NULL
    }
    
    inTrain <- rlply(k, createDataPartition(Y, p = p, list = FALSE))
    #inTrain <- rlply(k, sample(nr, round(p*nr)))
    outTrain <- lapply(inTrain, function(inTrain, total) total[-unique(inTrain)],
                       total = seq(nr))
    #i <- 1
    #method <- c("plsda")
    # loop through k bootstraps for stability metrics
    #require(tcltk)
    #pb <- tkProgressBar(title = "progress bar", min = 0,
    #                    max = k, width = 300)
    
    for(i in seq(k)){
      
      #setTkProgressBar(pb, i, label=paste( round(i/k*100, 0),
      #                                     "% done"))
      
      trainVars <- X[inTrain[[i]],, drop=F]
      trainVars.list[[i]] <- trainVars
      trainGroup <- Y[inTrain[[i]], drop=F]
      trainGroup.list[[i]] <- trainGroup
      
      trainData <- as.data.frame(trainVars)
      trainData$.classes <- trainGroup
      
      if(optimize == TRUE){
        if(optimize.resample == TRUE){
          # tune the methods
          tuned.methods <- optimize(trainVars = trainVars,
                                trainGroup = trainGroup,
                                method = method,
                                k.folds = k.folds,
                                repeats = repeats,
                                res = resolution,
                                grid = tuning.grid,
                                metric = metric,
                                savePerformanceMetrics = NULL,
                                allowParallel = allowParallel,
                                verbose = verbose,
                                theDots = theDots)
          
          # store the best tune parameters for each iteration
          names(tuned.methods$bestTune) = method
          resample.tunes[[i]] <- tuned.methods$bestTune
          
          if(i == 1){
            finalModel <- tuned.methods$finalModel
            finish.Model <- finalModel
          }else{
            finalModel <- tuned.methods$finalModel
            finish.Model <- append(finish.Model, tuned.methods$finalModel)
          }
          
          # Create empty list for features identified by each chosen algorithm
          features <- vector("list", length(method))
          names(features) <- tolower(method)
          
          for(j in seq(along = method)){
            ### Extract important features
            # pam requires a special mydata argument
            #mydata <- vector("list", length(method))
            if(method[j] == "pam"){
              #for(t in seq(along = method)){
                mydata <- list(x=t(trainVars.list[[i]]), y=factor(trainGroup.list[[i]]), geneid = as.character(colnames(trainVars.list[[i]])))
              #}
            }else{
              # svm requires training data for RFE
              #for(t in seq(along = method)){
                mydata <- trainVars.list[[i]]
              #}
            }
            
            features[[j]] <- extract.features(
              x = finalModel[j],
              dat = mydata,
              grp = trainGroup.list[[i]],
              # add in gbm best tune trees???
              bestTune = if(method[j] == "svm" | method[j] == "pam" | method[j] == "glmnet") tuned.methods$bestTune[[j]] else NULL,
              model.features = model.features, 
              method = method[j], 
              f = f, 
              comp.catch = tuned.methods$comp.catch)
            
            if(stability.metric %in% c("spearman", "canberra")){
              rownames(features[[j]]$features.selected) <- colnames(X)
            }
          }
          
          ### Re-fitting models to reduced features
          if(verbose){
            cat(paste("Refitting model iteration", i, "with selected features\n", sep = " "))
          }
          
          # subset only those features which were selected
          if(!is.null(f)){
            trainData.new <- lapply(features, FUN = function(x) trainData[,colnames(trainData) %in% c(t(x$features.selected), ".classes")])          
          }else{
            trainData.new <- lapply(features, FUN = function(x) trainData[,colnames(trainData) %in% c(rownames(x$features.selected), ".classes")])          
          }
          
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
                                          savePerformanceMetrics = FALSE,
                                          allowParallel = allowParallel,
                                          verbose = verbose,
                                          theDots = theDots)  
            }
          
          
          if(i == 1){
            finalModel.new <- sapply(tunedModel.new, FUN = function(x) x$finalModels)
            new.best.tunes <- sapply(tunedModel.new, FUN = function(x) x$bestTune)
            names(new.best.tunes) <- method
            final.features[[i]] <- sapply(features, FUN = function(x) x$features.selected)
            names(final.features[[i]]) <- method
          }else{
            tmp.model <- sapply(tunedModel.new, FUN = function(x) x$finalModels)
            tmp.tunes <- sapply(tunedModel.new, FUN = function(x) x$bestTune)
            names(tmp.tunes) <- method
            finalModel.new <- c(finalModel.new, tmp.model)
            new.best.tunes <- append(new.best.tunes, tmp.tunes)
            final.features[[i]] <- sapply(features, FUN = function(x) x$features.selected)
            names(final.features[[i]]) <- method
          }  
          
          # end of optimize.resample loop
        }else{
          if(i == 1){
            tuned.methods <- optimize(trainVars = trainVars,
                                  trainGroup = trainGroup,
                                  method = method,
                                  k.folds = k.folds,
                                  repeats = repeats,
                                  res = resolution,
                                  grid = tuning.grid,
                                  metric = metric,
                                  savePerformanceMetrics = NULL,
                                  verbose = verbose,
                                  allowParallel = allowParallel,
                                  theDots = theDots)
            finalModel <- tuned.methods$finalModel
            finish.Model <- tuned.methods$finalModel
          }else{
            # Fit remainder of resamples with initial best parameters
            tmp <- vector("list", length(method))
            names(tmp) <- method

            for(d in seq(along = method)){
              tmp[[d]] <- training(data = trainData,
                                   method = method[d],
                                   tuneValue = tuned.methods$bestTune[[d]],
                                   obsLevels = grp.levs,
                                   theDots = theDots)$fit
            }
            finalModel <- tmp
            finish.Model <- append(finish.Model, tmp)
          }
          
          # Create empty list for features identified by each chosen algorithm
          features <- vector("list", length(method))

          for(j in seq(along = method)){
            ### Extract important features
            # pam requires a special mydata argument
            if(method[j] == "pam"){
                mydata <- list(x=t(trainVars.list[[i]]), y=factor(trainGroup.list[[i]]), geneid = as.character(colnames(trainVars.list[[i]])))
            }else{
              # svm requires training data for RFE
              mydata <- trainVars.list[[i]]
            }
          
            features[[j]] <- extract.features(
              x = finalModel[j],
              dat = mydata,
              grp = trainGroup.list[[i]],
              # add in gbm best tune trees???
              bestTune = if(method[j] == "svm" | method[j] == "pam" | method[j] == "glmnet" | method[j] == "gbm") tuned.methods$bestTune[[j]] else NULL,
              model.features = model.features, 
              method = method[j], 
              f = f, 
              comp.catch = tuned.methods$comp.catch)
            
            if(stability.metric %in% c("spearman", "canberra")){
              rownames(features[[j]]$features.selected) <- colnames(X)
            }
          }
          
          ### Re-fitting models to reduced features
          if(verbose){
            cat(paste("Refitting model iteration", i, "with selected features\n", sep = " "))
          }
          
          # subset only those features which were selected
          if(!is.null(f) | model.features == TRUE){
            trainData.new <- lapply(features, FUN = function(x) trainData[,colnames(trainData) %in% c(t(x$features.selected), ".classes")])          
          }else{
            trainData.new <- lapply(features, FUN = function(x) trainData[,colnames(trainData) %in% c(rownames(x$features.selected), ".classes")])          
          }
          
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
                                          savePerformanceMetrics = FALSE,
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
                                   tuneValue = as.data.frame(t(unlist(tunedModel.new[[d]]$bestTune))),
                                   obsLevels = grp.levs,
                                   theDots = theDots)$fit
            }
          }

          if(i == 1){
            finalModel.new <- sapply(tunedModel.new, FUN = function(x) x$finalModels)
            new.best.tunes <- sapply(tunedModel.new, FUN = function(x) x$bestTune)
            final.features[[i]] <- sapply(features, FUN = function(x) x$features.selected)
            #names(final.features[[i]]) <- method
          }else{
            tmp.model <- lapply(tmp, FUN = function(x) x)
            finalModel.new <- c(finalModel.new, tmp.model)
            final.features[[i]] <- sapply(features, FUN = function(x) x$features.selected)
            #names(final.features[[i]]) <- method
          }  
          
        } # end of single optimized loop
        # end of optimize loop 
      }else{        
        #theDots <- list(ncomp = 3, mtry = 10)
        names(theDots) <- paste(".", names(theDots), sep="")
        
        # sequester appropriate parameters to fit models
        args.seq <- sequester(theDots, method)
        
        # remove arguments used from theDots - also remove '.' from each
        names(theDots) <- sub(".", "", names(theDots))
        moreDots <- theDots[!names(theDots) %in% args.seq$pnames]
        if(length(moreDots) == 0){
          moreDots <- NULL
        }
                
        tmp <- vector("list", length(method))
        for(q in seq(along = method)){
          tmp[[q]] <- training(data = trainData,
                                      method = method[q],
                                      tuneValue = args.seq$parameters[[q]],
                                      obsLevels = grp.levs,
                                      theDots = moreDots)$fit  
        }
        
        finalModel <- tmp
        if(i == 1){
          finish.Model <- finalModel
        }else{
          finish.Model <- append(finish.Model, tmp)
        }
        
        # Create empty list for features identified by each chosen algorithm
        features <- vector("list", length(method))
        names(features) <- tolower(method)
        
        for(j in seq(along = method)){
          ### Extract important features
          # pam requires a special mydata argument
          mydata <- vector("list", length(method))
          if(method[j] == "pam"){
            for(t in seq(along = method)){
              mydata[[t]] <- list(x=t(trainVars.list[[i]]), y=factor(trainGroup.list[[i]]), geneid = as.character(colnames(trainVars.list[[i]])))
            }
          }else{
            # svm requires training data for RFE
            for(t in seq(along = method)){
              mydata[[t]] <- trainVars.list[[i]]
            }
          }
          
          features[[j]] <- extract.features(
            x = finalModel[j],
            dat = mydata[[j]],
            grp = trainGroup.list[[j]],
            # add in gbm best tune trees???
            bestTune = if(method[j] == "svm" | method[j] == "pam" | method[j] == "glmnet") args.seq$parameters[[j]] else NULL,
            model.features = model.features, 
            method = method[j], 
            f = f, 
            comp.catch = tuned.methods$comp.catch)
          
          if(stability.metric %in% c("spearman", "canberra")){
            rownames(features[[j]]$features.selected) <- colnames(X)
          }
        }
        
        ### Re-fitting models to reduced features
        if(verbose){
          cat(paste("Iteration ", i, " Refitting models with selected features\n", sep = ""))
        }
        
        # subset only those features which were selected
        if(!is.null(f)){
          trainData.new <- lapply(features, FUN = function(x) trainData[,colnames(trainData) %in% c(t(x$features.selected), ".classes")])          
        }else{
          trainData.new <- lapply(features, FUN = function(x) trainData[,colnames(trainData) %in% c(rownames(x$features.selected), ".classes")])          
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
          finalModel.new <- lapply(tunedModel.new, FUN = function(x) x)
          names(finalModel.new) <- method
          final.features[[i]] <- sapply(features, FUN = function(x) x$features.selected)
          names(final.features[[i]]) <- method
        }else{
          tmp.model <- lapply(tunedModel.new, FUN = function(x) x)
          names(tmp.model) <- method
          finalModel.new <- c(finalModel.new, tmp.model)
          final.features[[i]] <- sapply(features, FUN = function(x) x$features.selected)
          names(final.features[[i]]) <- method
        }  
      } # end of non-optimized sequence
    } # end of stability loop
        
    #close(pb)
    
    ### Performance Metrics of Reduced Models
    #if(verbose){
      cat("Calculating Model Performance Statistics\n")
    #}
    
    
    final.metrics <- prediction.metrics(finalModel = finalModel.new,
                                        method = method,
                                        raw.data = raw.data,
                                        inTrain = inTrain,
                                        outTrain = outTrain,
                                        features = final.features,
                                        bestTune = if(optimize) new.best.tunes else args.seq$parameters,
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
          tmp <- subset(x, rownames(x) == method[h])
          performance[[h]] <- c(colMeans(tmp, na.rm = TRUE), apply(x, 2, sd, na.rm = TRUE))
          names(performance[[h]])[-(1:ncol(tmp))] <- paste(names(performance[[h]])[-(1:ncol(tmp))], "SD", sep = " ")
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
    #c <- 3
    if(!model.features){
      for(c in seq(along = method)){
        #met <- method[c]
        #if(met == "svm" | met == "glmnet"){
        #  results.stability[[c]] <- as.data.frame(sapply(final.features, FUN = function(x) x))
        #}else{
        results.stability[[c]] <- as.data.frame(sapply(final.features, FUN = function(x) x[[c]]))
        #}
        if(is.null(f)){
          rownames(results.stability[[c]]) <- colnames(X)
        }
      }
    }else{
      for(c in seq(along = method)){
        tmp <- lapply(final.features, FUN = function(x) as.data.frame(t(data.frame(x[[c]]))))
        results.stability[[c]] <- t(rbind.fill(tmp))
        rownames(results.stability[[c]]) <- NULL
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
    
    # add stability metrics to features selected
    for(i in seq(along = method)){
      results.stability[[i]] <- list(metric = stability.metric, features = results.stability[[i]], stability = stability[[i]])
    }
    
    ### Desired Output 
    ## specifications
    specs = list(total.samples=nr, 
                 number.features=nc, 
                 number.groups=num.group, 
                 group.levels=grp.levs, 
                 number.observations.group=num.obs.group)
    
    ## add remainder of data
    overall <- list(methods = method,                     # algorithms run
                    performance = performance,            # performance metrics of each algorithm
                    RPT = rpt,                             # robustness-performance trade-off
                    features = results.stability,         # list of each algorithms results
                    stability.models = stability.models,  # stability amongst algorithms
                    original.best.tunes = resample.tunes, # if optimize.resample returns the best tunes for each iteration
                    final.best.tunes = if(optimize.resample) all.model.perfs else NULL,   # if optimize.resample, provide parameter with performance statistics
                    specs = specs                         # general specs of the input data
    )
    return(overall)
  }

    