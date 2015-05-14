
#' @title Fit Models without Feature Selection
#' @description Applies models to high-dimensional data for classification.
#' @param X A scaled matrix or dataframe containing numeric values of 
#' each feature
#' @param Y A factor vector containing group membership of samples
#' @param method A vector listing models to be fit.
#' Available options are \code{"plsda"} (Partial Least Squares Discriminant 
#' Analysis), \code{"rf"} (Random Forest), \code{"gbm"} (Gradient Boosting 
#' Machine), \code{"svm"} (Support Vector Machines), \code{"glmnet"} 
#' (Elastic-net Generalized Linear Model), and \code{"pam"} (Prediction 
#' Analysis of Microarrays)
#' @param p Percent of data to by 'trained'
#' @param optimize Logical argument determining if each model should be 
#' optimized. Default \code{"optimize = TRUE"}
#' @param tuning.grid Optional list of grids containing parameters to optimize 
#' for each algorithm. Default \code{"tuning.grid = NULL"} lets function 
#' create grid determined by \code{"res"}
#' @param k.folds Number of folds generated during cross-validation.  
#' Default \code{"k.folds = 10"}
#' @param repeats Number of times cross-validation repeated.  
#' Default \code{"repeats = 3"}
#' @param resolution Resolution of model optimization grid.  
#' Default \code{"resolution = 3"}
#' @param metric Criteria for model optimization.  
#' Available options are \code{"Accuracy"} (Predication Accuracy),
#' \code{"Kappa"} (Kappa Statistic), and \code{"AUC-ROC"} (Area Under the 
#' Curve - Receiver Operator Curve)
#' @param allowParallel Logical argument dictating if parallel processing 
#' is allowed via foreach package.
#' Default \code{allowParallel = FALSE}
#' @param verbose Logical argument if should output progress
#' @param ... Extra arguments that the user would like to apply to the models
#'
#' @return \item{Methods}{Vector of models fit to data}
#' @return \item{performance}{Performance metrics of each model and 
#' bootstrap iteration}
#' @return \item{specs}{List with the
#' following elements:}
#' @return \itemize{
#'  \item{total.samples: Number of samples in original dataset}
#'  \item{number.features: Number of features in orginal dataset}
#'  \item{number.groups: Number of groups}
#'  \item{group.levels: The specific levels of the groups}
#'  \item{number.observations.group: Number of observations in each group}}
#' @author Charles Determan Jr
#' @example inst/examples/fit.only.model.R
#' @import DiscriMiner
#' @import randomForest
#' @import e1071
#' @import gbm
#' @import pamr
#' @import glmnet
#' @export

fit.only.model <- 
    function(X,
             Y,
             method,
             p = 0.9,
             optimize = TRUE,
             tuning.grid = NULL,
             k.folds = if(optimize) 10 else NULL,
             repeats = if(optimize) 3 else NULL,
             resolution = if(optimize) 3 else NULL,
             metric = "Accuracy",
             allowParallel = FALSE,
             verbose = 'none',
             ...
    )
{    
        #### Filter Methods???
        ## Bhattacharyya Distance
        #??bhattacharyya
        #require(fpc)
        #bhattacharyya.dist
        ## Relative Entropy
        
        assert_is_character(verbose)
        
        verify(x = X, y = Y, method = method, na.rm = FALSE)
        
        if (is.null(colnames(X))){
            colnames(X) = paste(rep("X",ncol(X)),seq_len(ncol(X)), sep='') 
        }
        if (is.null(rownames(X))) rownames(X) = 1:nrow(X)
        
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
        
        inTrain <- caret::createDataPartition(Y, p = p, list = FALSE)
        outTrain <- seq(nr)[-unique(inTrain)]
        trainVars <- X[inTrain,, drop=FALSE]
        trainGroup <- Y[inTrain, drop=FALSE]
        trainData <- as.data.frame(trainVars)
        trainData$.classes <- trainGroup
        
        # must add repeat loops (same as k in fs.stability for StDev)
        
        if(optimize == TRUE){
            tuned.methods <- 
                optimize.model(
                    trainVars = trainVars,
                    trainGroup = trainGroup,
                    method = method,
                    k.folds = k.folds,
                    repeats = repeats,
                    res = resolution,
                    grid = tuning.grid,
                    metric = metric,
                    #savePerformanceMetrics = NULL,
                    verbose = verbose,
                    allowParallel = allowParallel,
                    theDots = theDots)
            finalModel <- tuned.methods$finalModel
            best.tunes <- tuned.methods$bestTune
            names(best.tunes) <- method
            
            # end of single optimized loop
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
            
            finish.Model <- tmp
            
        } # end of non-optimized sequence
        
        ### Performance Metrics
        cat("Calculating Model Performance Statistics\n")    
        
        final.metrics <- 
            prediction.metrics(
                finalModel = finalModel,
                method = method,
                raw.data = raw.data,
                inTrain = inTrain,
                outTrain = outTrain,
                features = NULL,
                bestTune = if(optimize) best.tunes else args.seq$parameters,
                grp.levs = grp.levs,
                stability.metric = NULL)
        
        ### Extract Performance Metrics
        if(optimize == TRUE){
            colnames(final.metrics) <- gsub("^\\.", "", colnames(final.metrics))  
            performance <- vector("list", length(method))
            names(performance) <- method
            x <- final.metrics[,!grepl("^cell", 
                                       colnames(final.metrics)),
                               drop = FALSE]
            for(h in seq(along = method)){
                tmp <- subset(x, rownames(x) == method[h])
                performance[[h]] <- c(colMeans(tmp, na.rm = TRUE), 
                                      apply(tmp, 2, sd, na.rm = TRUE))
                names(performance[[h]])[-(1:ncol(tmp))] <- 
                    paste(names(performance[[h]])[-(1:ncol(tmp))], 
                          "SD", 
                          sep = " ")
                performance[[h]] <- 
                    do.call(cbind, 
                            c(as.vector(best.tunes[[h]]), performance[[h]]))
                colnames(performance[[h]]) <- 
                    gsub("^\\.", "", colnames(performance[[h]]))
                rownames(performance[[h]]) <- 1
            }
        }else{
            colnames(final.metrics) <- 
                gsub("^\\.", "", colnames(final.metrics))  
            performance <- vector("list", length(method))
            names(performance) <- method
            for(h in seq(along = method)){
                x <- final.metrics[,!grepl("^cell", colnames(final.metrics)),
                                   drop = FALSE]
                tmp <- subset(x, rownames(x) == method[h])
                performance[[h]] <- c(colMeans(x, na.rm = TRUE), 
                                      apply(tmp, 2, sd, na.rm = TRUE))
                names(performance[[h]])[-(1:ncol(tmp))] <- 
                    paste(names(performance[[h]])[-(1:ncol(tmp))], 
                          "SD", sep = " ")
                performance[[h]] <- 
                    do.call(cbind, 
                            c(as.vector(args.seq$parameters[[h]]), 
                              performance[[h]]
                              )
                            )
                colnames(performance[[h]]) <- 
                    gsub("^\\.", "", colnames(performance[[h]]))
                rownames(performance[[h]]) <- 1
            }
        }
        
        
        ### Desired Output 
        ## specifications
        specs = list(total.samples=nr, 
                     number.features=nc, 
                     number.groups=num.group, 
                     group.levels=grp.levs, 
                     number.observations.group=num.obs.group)
        
        ## add remainder of data
        overall <- list(methods = method,
                        performance = performance,
                        specs = specs
        )
        return(overall)
    }
