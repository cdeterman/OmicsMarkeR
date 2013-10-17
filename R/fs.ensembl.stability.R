
##################################
### Ensemble Feature Selection ###
##################################

fs.ensembl.stability <- 
  function(variables,                          # scaled matrix or dataframe of explanatory variables
           groups,                             # vector or factor with group memberships
           method,                             # "PLSDA", "glmnet","SVM", "RF", "GBM","PAM"
           k = 10,                             # number of subsamples
           p = 0.9,                            # percentage data subsampled
           f = ceiling(ncol(variables)/10),    # number of top features subset
           bags=40,                            # number of bags for bootstrap aggregation (i.e. bagging)
           aggregation.metric = "CLA",         # Method to aggregation ensemble results
           stability.metric = "jaccard",       # stability metric
           optimize = TRUE,
           optimize.resample = FALSE,
           tuning.grid = NULL,
           k.folds = if(optimize) 10 else NULL,
           repeats = if(optimize) 3 else NULL,
           resolution = if(optimize) 3 else NULL,
           metric = "Accuracy",
           model.features = FALSE,
           verbose = FALSE,
           ...
           )
    {    
    ### load all libraries
    #if("PLSDA" %in% method) # require(DiscriMiner) removed until package updated
    #if("PAM" %in% method) require(pamr)
    #if("glmnet" %in% method) require(glmnet)
    #if("RF" %in% method) require(randomForest)
    #if("GBM" %in% method) require(gbm)
    #if("SVM" %in% method) require(e1071) # or require(kernlab)
    
    verify_data <- verify(x = variables, y = groups, method = method, f = f, stability.metric = stability.metric, model.features = model.features, na.rm = FALSE)
    #verify_data <- my_verify(variables, groups, na.rm = FALSE)
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
    
    inTrain <- rlply(k, sample(nr, round(p*nr)))
    outTrain <- lapply(inTrain, function(inTrain, total) total[-unique(inTrain)],
                       total = seq(nr))
    #i <- 1
    #optimize.resample = FALSE
    
    
    ### Stability Loop
    for (i in 1:k){
      # random sample of data
      #inTrain <- sample(nr, round(p*(nr)))
      
      trainX <- X[inTrain[[i]],, drop=F]
      trainY <- Y[inTrain[[i]], drop=F]
      trainData <- as.data.frame(trainVars)
      trainData$.classes <- trainY
      
      ## Bagging loop
      results.bagging <- bagging.wrapper(X = trainX, 
                                         Y = trainY, 
                                         method = method,
                                         bags = bags, 
                                         f = f, 
                                         aggregation.metric = aggregation.metric,
                                         #inTrain = inTrain[[i]],
                                         k.folds = k.folds,
                                         repeats = repeats,
                                         res = res,
                                         tuning.grid = tuning.grid,
                                         optimize = optimize,
                                         optimize.resample = optimize.resample,
                                         metric = metric,
                                         model.features = model.features,
                                         verbose = verbose,
                                         theDots = theDots)
      
      # store the best tune parameters for each iteration
      if(!is.null(bagged.tunes)){
        bagged.tunes[[i]] <- results.bagging$bestTunes
      }
      
      # Store the features selected for stability analysis
      features[[i]] <- results.bagging$results$ensemble.results
      
      ### Re-fitting models to reduced features
      # subset only those features which were selected
      trainData.new <- lapply(features[[i]], FUN = function(x) trainData[,colnames(trainData) %in% c(as.vector(names(x)), ".classes")])
      
      if(optimize == TRUE){
        if(optimize.resample == TRUE){
          tunedModel.new <- vector("list", length(method))
          for(m in seq(along = method)){
            tunedModel.new[[m]] <- tune(trainVars = trainData.new[[m]][,!colnames(trainData.new[[m]]) %in% c(".classes")],
                                        trainGroup = trainData.new[[m]]$.classes,
                                        method = method[m],
                                        k.folds = k.folds,
                                        repeats = repeats,
                                        res = resolution,
                                        grid = tuning.grid,
                                        metric = metric,
                                        savePerformanceMetrics = FALSE,
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
            finalModel.new <- append(finalModel.new, tmp.model)
            new.best.tunes <- append(new.best.tunes, tmp.tunes)
          }  
          # end of full optimize loop
        }else{
          if(i == 1){
            #m <- 1
            tunedModel.new <- vector("list", length(method))
            for(m in seq(along = method)){
              tunedModel.new[[m]] <- tune(trainVars = trainData.new[[m]][,!colnames(trainData.new[[m]]) %in% c(".classes")],
                                          trainGroup = trainData.new[[m]]$.classes,
                                          method = method[m],
                                          k.folds = k.folds,
                                          repeats = repeats,
                                          res = resolution,
                                          grid = tuning.grid,
                                          metric = metric,
                                          savePerformanceMetrics = FALSE,
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
              #tunedModel.new[d]
            #length(tunedModel.new)
            #  data.frame(t(unlist(tunedModel.new[[1]]$bestTune)))$.ncomp
            #  str(tunedModel.new)
            }
          }
          
          if(i == 1){
            finalModel.new <- sapply(tunedModel.new, FUN = function(x) x$finalModels)
            new.best.tunes <- sapply(tunedModel.new, FUN = function(x) x$bestTune)
          }else{
            tmp.model <- sapply(tmp, FUN = function(x) x)
            finalModel.new <- append(finalModel.new, tmp.model)
          }  
        } # end of single optimize loop
        # end of optimize loop
      }else{
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
          tmp.model <- sapply(tunedModel.new, FUN = function(x) x)
          finalModel.new <- append(finalModel.new, tmp.model)
        }
      } # end of non-optimized loop
    } # end stability loop
    
    ### Performance Metrics of Reduced Models
    final.metrics <- prediction.metrics(finalModel = finalModel.new,
                                        method = method,
                                        raw.data = raw.data,
                                        inTrain = inTrain,
                                        outTrain = outTrain,
                                        bestTune = new.best.tunes,
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
    stability <- lapply(results.stability, pairwise.stability, stability.metric = stability.metric, k = k, nc= nc)
    
    # stability across algorithms (i.e. 'function perturbation' ensemble analysis)
    if(length(method) > 1){
      stability.models <- pairwise.model.stability(ranked_features = results.stability,
                                                   stability.metric = stability.metric,
                                                   m = length(method),
                                                   k = k)
    }else{
      stability.models <- NULL
    }
    
    # add stability metric to each respective algorithm
    for(i in 1:length(method)){
      results.stability[[i]] <- list(metric = stability.metric, results.stability = results.stability[[i]], stability = stability[[i]])
    }
    
    specs = list(total.samples=nr, 
                 number.features=nc, 
                 number.groups=num.group, 
                 group.levels=grp.levs, 
                 number.observations.group=num.obs.group)
    
    ## add remainder of data
    overall <- structure(list(methods = method,                      # algorithms run
                              performance = performance,
                              features = results.stability,          # list of each algorithms results
                              stability.models = stability.models,   # stability amongst algorithms
                              all.tunes = bagged.tunes,              # if optimize.resample returns the best tunes for each iteration
                              final.best.tunes = if(optimize.resample) all.model.perfs else NULL,   # if optimize.resample, provide parameter with performance statistics
                              specs = specs                          # general specs of the input data
    ),
                         class = "ensemble.stability")
    return(overall)
}
