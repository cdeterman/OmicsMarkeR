
#' @title Prediction Metric Calculations
#' @description Performance evaluation of all fitted models.  This function 
#' concisely provides model performance metrics, including confusion matrix 
#' and ROC.
#' @param finalModel List of fitted models
#' @param method Vector of strings dictating the models that were fit
#' @param raw.data Original dataset prior to any training subset
#' @param inTrain List of training indicies for each feature selection run
#' @param outTrain List of testing data indicies for each feature selection run
#' @param features List of selected features for each model
#' @param bestTune List of parameters that have been optimized for the each 
#' respective model
#' @param grp.levs Vector of group levels
#' @param stability.metric A character object specifying the stability metric
#' @return Returns a dataframe consisting of each feature selection runs 
#' evaluated Accuracy, Kappa, ROC.AUC, Sensitivity, Specificity, Positive 
#' Predictive Value, and Negative Predictive Value.
#' @seealso \code{\link{performance.stats}}, 
#' \code{\link{perf.calc}} caret function \link[caret]{confusionMatrix}
# ' @export

prediction.metrics <- 
    function(finalModel,
             method,
             raw.data,
             inTrain,
             outTrain,
             features,
             bestTune,
             grp.levs,
             stability.metric)
    {
        raw.data.vars <- raw.data[,!colnames(raw.data) %in% c(".classes")]
        raw.data.grps <- raw.data[,colnames(raw.data) %in% c(".classes")]
        
        if(class(inTrain) == "list" & class(outTrain) == "list"){
            inTrain.list <- rep(inTrain, length(method))
            outTrain.list <- rep(outTrain, length(method))
        }else{
            inTrain.list <- rep(list(inTrain), length(finalModel))
            outTrain.list <- rep(list(outTrain), length(finalModel))
        }
        
        # a check for when results come from optimize.resample = FALSE
        if(length(bestTune) != length(finalModel)){
            tmp.mult <- length(finalModel)/length(bestTune)
            bestTune <- rep(bestTune, tmp.mult)
            names(bestTune) <- names(finalModel)
        }
        
        # sort the list elements so applied in the proper order
        method.names <- 
            unlist(
                lapply(method, 
                       FUN = function(x){
                           c(rep(x, 
                                 length(bestTune)/length(method)))
                       }
                )
            )
        bestTune <- bestTune[match(method.names,names(bestTune))]
        finalModel <- finalModel[match(method.names, names(finalModel))]    
        
        # check features
        if(is.null(features)){
            features <- vector("list", length(finalModel))
            for(f in seq(length(finalModel))){
                features[[f]] <- colnames(raw.data.vars)
            }
        }
        
        features <- features[match(method.names, names(features))]
        
        predicted <- vector("list", length(finalModel))
        names(predicted) <- names(finalModel)
        
        for(e in seq(along = finalModel)){
            new.dat <- 
                switch(names(finalModel[e]),
                       svm = {
                           if(stability.metric %in% c("spearman", "canberra")){
                               raw.data.vars[outTrain.list[[e]],, 
                                             drop = FALSE]
                           }else{
                               raw.data.vars[outTrain.list[[e]],
                                             (names(raw.data.vars) 
                                              %in% features[[e]]), 
                                             drop = FALSE]
                           }
                           },
                       
                       glmnet= {
                           if(stability.metric %in% c("spearman", "canberra")){
                               raw.data.vars[outTrain.list[[e]],, 
                                             drop = FALSE]
                           }else{
                               raw.data.vars[outTrain.list[[e]],
                                             (names(raw.data.vars) 
                                              %in% features[[e]]), 
                                             drop = FALSE]
                           }
                       },
                       
                       pam = {
                           if(stability.metric %in% c("spearman", "canberra")){
                               raw.data.vars[outTrain.list[[e]],, 
                                             drop = FALSE]
                           }else{
                               raw.data.vars[outTrain.list[[e]],
                                             (names(raw.data.vars) 
                                              %in% features[[e]]), 
                                             drop = FALSE]
                           }
                       },
                       
                       plsda =, gbm =, rf = {raw.data.vars[outTrain.list[[e]],
                                                           ,drop = FALSE]},
                )
            
            predicted[[e]] <- predicting(method = names(finalModel)[e],
                                         modelFit = finalModel[[e]],
                                         orig.data = raw.data,
                                         indicies = inTrain.list[[e]],
                                         newdata = new.dat,
                                         param = bestTune[[e]]
            )
        }
        
        ## collate the predicitons across all the models
        for(g in seq(along = finalModel)){
            predicted[[g]] <- 
                factor(as.character(unlist(predicted[[g]])), levels = grp.levs)
            predicted[[g]] <- 
                data.frame(pred = predicted[[g]], 
                           obs = raw.data.grps[outTrain.list[[g]]], 
                           stringsAsFactors = FALSE)
        }
        
        method.vector <- rep(method, each = length(finalModel)/length(method))
        perf.metrics <- mapply(predicted,
                               FUN = function(x, y) perf.calc(x, 
                                                              lev = grp.levs, 
                                                              model = y),
                               y = method.vector, SIMPLIFY = FALSE)                 
        
        cells <- lapply(predicted,
                        function(x) flatTable(x$pred, x$obs))
        for(ind in seq(along = cells)){
            perf.metrics[[ind]] <- c(perf.metrics[[ind]], cells[[ind]])
        } 
        
        final.metrics <- do.call("rbind", perf.metrics)
    }
