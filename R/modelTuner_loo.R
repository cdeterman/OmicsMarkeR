
# declare global variables (i.e. the foreach iterators)
globalVariables(c('algo', 'parms', 'iter'))

#' @title Model Tuner for Leave-One-Out Cross-Validation
#' @description Optimizes each model via LOO CV based upon the parameters 
#' provided either by the internal \code{\link{denovo.grid}}
#' function or by the user.
#' @param trainData Data used to fit the model
#' @param guide Output from \code{\link{tune.instructions}}.  Facilitates the 
#' optimization by avoiding redundant model fitting.
#' @param method Vector of strins listing models to be fit
#' @param inTrain Indicies for cross-validated training models
#' @param outTrain Indicies for cross-validated testing models
#' @param lev Group levels
#' @param savePredictions Logical argument dictating if should save the 
#' prediction data.  Default \code{savePredictions = FALSE}
#' @param allowParallel Logical argument dictating if parallel processing is 
#' allowed via foreach package
#' @param verbose Logical argument dictating if should print progress
#' @param theDots List of additional arguments provided in the initial 
#' classification and features selection function
#' @return Returns list of fitted models
#' @import DiscriMiner
#' @import randomForest
#' @import plyr
#' @importFrom caret progress
#' @import e1071
#' @import gbm
#' @import pamr
#' @import glmnet
#' @import foreach
#' @import caTools
# ' @export
#' @author Charles E. Determan Jr.

modelTuner_loo <- function(trainData,
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
    # Set parallel option  
    `%op%` <- if(allowParallel){
        `%dopar%` 
    }else{
        `%do%`
    }  
    
    # set number of iterations to number of resample lists 
    # (i.e. number of samples)
    # iter = seq(along = inTrain)
    
    # This is a triple loop
    # first loop cycles through the methods chosen
    # second loop cycles through of the LOO cross validations
    # third loop cycles through parameters of each cross validation cycle
    
    tmp.list <- foreach(algo = seq(along = method),
                        .verbose = FALSE, 
                        .packages = c("OmicsMarkeR","foreach","plyr",
                                      "DiscriMiner","randomForest","e1071",
                                      "gbm","pamr","glmnet","caTools"),
                        .export = c("training", "round.multiple", 
                                    "predicting", "expandParameters",
                                    "flatTable"),
                        .errorhandling = "stop") %:%
        foreach(iter = seq(along = inTrain), 
                .combine = "c", 
                .verbose = FALSE, 
                .packages = c("OmicsMarkeR","foreach","plyr","DiscriMiner",
                              "randomForest","e1071","gbm","pamr","glmnet",
                              "caTools"),
                .export = c("training", "round.multiple", "predicting", 
                            "expandParameters", "flatTable"),
                .errorhandling = "stop") %:%
        foreach(parms = seq(nrow(guide[[algo]]$loop)),
                .combine = "c", 
                .verbose = FALSE, 
                .packages = c("OmicsMarkeR","foreach","plyr","DiscriMiner",
                              "randomForest","e1071","gbm","pamr","glmnet",
                              "caTools"),
                .export = c("training", "round.multiple", "predicting", 
                            "expandParameters", "flatTable"),
                .errorhandling = "stop") %op%
{      
    ## Removes '.' from start of each parameter
    ## create 'printed' for verbose printing
    printed <- format(guide[[algo]]$loop, digits = 4)
    colnames(printed) <- gsub("^\\.", "", colnames(printed))
    
    # library(caret)
    if(verbose) progress(printed[parms,,drop = FALSE],
                         names(inTrain), iter)
    
    #if(testing) cat("pre-model\n")
    outIndex <- outTrain[[iter]]
    
    # create models
    mod <- training(data = trainData[complete.cases(
        trainData[inTrain[[iter]],,drop = FALSE]),,drop = FALSE],
        method = method[algo],
        tuneValue = guide[[algo]]$loop[parms,,drop = FALSE],
        obsLevels = lev,
        theDots = theDots)
    
    predicted <- 
        predicting(method = method[algo],
                   modelFit = mod$fit,
                   orig.data = trainData,
                   indicies = inTrain[[iter]],
                   newdata = trainData[outIndex, 
                                       !(names(trainData) %in% ".classes"), 
                                       drop = FALSE],
                   param = guide[[algo]]$seqParam[[parms]])
    
    ##################################
    
    # If the model was built with parameters that 'submodels' 
    # can be extracted this section will combine them together
    if(!is.null(guide[[algo]]$seq))
    {
        
        ## merge the fixed and seq parameter values together
        allParam <- expandParameters(guide[[algo]]$loop[parms,,drop = FALSE], 
                                     guide[[algo]]$seqParam[[parms]])
        
        predicted <- lapply(predicted,
                            function(x, y, lv)
                            {
                                if(!is.factor(x) & is.character(x)){
                                    x <- factor(as.character(x), levels = lv)
                                } 
                                data.frame(pred = x, 
                                           obs = y, 
                                           stringsAsFactors = FALSE)
                            },
                            y = trainData$.classes[outIndex],
                            lv = lev)
        
        predicted <- do.call("rbind", predicted)
        predicted <- cbind(predicted, allParam)
        
    } else {
        
        # for models without retaining 'lower' parameters
        if(is.factor(trainData$.classes)){
            predicted <- factor(as.character(predicted),
                                levels = lev)
        } 
        predicted <-  data.frame(pred = predicted,
                                 obs = trainData$.classes[outIndex],
                                 stringsAsFactors = FALSE)
        predicted <- cbind(predicted, guide[[algo]]$loop[parms,,drop = FALSE])
    }
    predicted$sampleIndex <- names(inTrain)[iter]
    
    # Print Progress
    if(verbose) progress(printed[parms,,drop=FALSE],
                         names(inTrain), iter, FALSE)
    list(tunes=predicted)
}

####################################
###### Tuning loops Complete #######
####################################
if(length(method) > 1){
    names(tmp.list) <- method
    
    ## plyr:::rbind.fill - binds list of dataframes together
    tunes <- lapply(tmp.list, FUN = function(x) 
        rbind.fill(x[names(x) == "tunes"]))
    
    ## remove '.' from each name
    new.names <- lapply(tunes, FUN = function(x) 
        gsub("^\\.", "", names(x)))
    tunes <- mapply(tunes, FUN = function(x,y) {
        names(x) <- y; return(x)}, 
        y = new.names, SIMPLIFY = FALSE)
    
    for(i in length(tunes)){
        if(any(
            !complete.cases(
                tunes[[i]][,!grepl("^cell|sampleIndex", names(tunes[[i]])),
                           drop = FALSE]
            )
        ))
        {
            warning(paste("There were missing values in resampled 
                          performance measures in", 
                          names(tunes[i]), sep = " "))
        } 
        }
    
    # Retaining the parameters
    par_levs <- lapply(tunes, 
                       FUN = function(x) 
                           unique(
                               sapply(x[,!colnames(x) %in% 
                                               c("pred", "obs", "sampleIndex"), 
                                           drop = FALSE], as.factor)
                               )
                       )
    
    # split each 'method' list into multiple lists
    split.results <- vector("list", length(method))
    names(split.results) <- method
    
    for(im in 1:length(method)){
        lookup <- params(method[im])[[method[im]]]$parameter
        #lookup <- params("gbm")$gbm$parameter
        if(length(lookup) > 1){
            filter <- paste(lookup,sep= ",")
        }else{
            filter <- as.character(lookup)
        }
        
        # Get performance metrics
        split.results[[im]] <- do.call("rbind", lapply(
            split(tunes[[im]], tunes[[im]][,c(filter)]),
            FUN = perf.calc,
            lev = lev))  
        # Make sure parameters are numeric for subsequent sorting
        # Bind parameters to results
        # Rename rows to make pretty
        split.results[[im]] <- 
            cbind(sapply(
                as.data.frame(par_levs[[im]]), 
                FUN = function(x) as.numeric(as.character(x))), 
                split.results[[im]])
        rownames(split.results[[im]]) <- seq(nrow(split.results[[im]]))
    }
    
    print("Model Tuning Complete")
    
    out <- vector("list", length(method))
    names(out) <- method
    for(i in seq(along = method)){
        out[[i]] <- list(performance = split.results[[i]], tunes = tunes[[i]])
    }
    }else{
        
        tmp.list <- unlist(tmp.list, recursive = FALSE)
        ## plyr:::rbind.fill - binds list of dataframes together
        tunes <- rbind.fill(tmp.list[names(tmp.list) == "tunes"])
        
        ## remove '.' from each name
        names(tunes) <- gsub("^\\.", "", names(tunes))  
        if(any(!complete.cases(
            tunes[,!grepl("^cell|sampleIndex", 
                          colnames(tunes)),
                  drop = FALSE])))
        {
            warning("There were missing values in resampled 
                    performance measures.")
        }
        
        par_levs <- unique(
            sapply(
                tunes[,!colnames(tunes) %in% c("pred", "obs", "sampleIndex"), 
                      drop = FALSE], as.factor))
        
        lookup <- params(method)[[method]]$parameter
        if(length(lookup) > 1){
            filter <- paste(lookup, sep=",")
        }else{
            filter <- as.character(lookup)
        }
        
        # Get performance metrics
        # make sure numeric for subsequent sorting
        # bind metrics and parameters together
        split.results <- do.call("rbind", 
                                 lapply(
                                     split(tunes, tunes[,c(filter)]),
                                     FUN = perf.calc,
                                     lev = lev)
                                 )
        split.results <- cbind(
            sapply(as.data.frame(par_levs), FUN = function(x){
                as.numeric(as.character(x))}), 
            split.results)
        rownames(split.results) <- seq(nrow(split.results))
        
        print(paste(method, "complete"))
        out <- list(performance = split.results, tunes = tunes)
    }

out
}