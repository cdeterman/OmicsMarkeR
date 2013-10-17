
### Prediction Metric Calculations

prediction.metrics <- 
  function(finalModel,
           method,
           raw.data,
           inTrain,
           outTrain,
           bestTune,    #tuned.methods$bestTune
           grp.levs)
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
    }
    
    # if list element contains parameters of model, name the element of that model
    test.param <- lapply(params(method), FUN = function(x) paste(".", as.character(x$parameter), sep=""))
    for(i in seq(along = bestTune)){
      model.index <- which(test.param == names(bestTune[[i]]))
      names(bestTune)[[i]] <- names(model.index)
    }    
    # sort the list elements so applied in the proper order
    method.names <- unlist(lapply(method, FUN = function(x) c(rep(x, length(bestTune)/length(method)))))
    bestTune <- bestTune[order(names(bestTune), levels = method.names)]
    finalModel <- finalModel[order(names(finalModel), levels = method.names)]
    
    predicted <- vector("list", length(finalModel))
    names(predicted) <- names(finalModel)
    
    for(e in seq(along = finalModel)){
      predicted[[e]] <- predicting(method = names(finalModel)[e],
                                   modelFit = finalModel[[e]],
                                   orig.data = raw.data,
                                   indicies = inTrain.list[[e]],
                                   newdata = raw.data.vars[outTrain.list[[e]],, drop = FALSE],
                                   param = bestTune[[e]]
      )
    }
    
    
    ## collate the predicitons across all the models
    for(g in seq(along = finalModel)){
      predicted[[g]] <- factor(as.character(unlist(predicted[[g]])), levels = grp.levs)
      predicted[[g]] <- data.frame(pred = predicted[[g]], obs = raw.data.grps[outTrain.list[[g]]], stringsAsFactors = FALSE)
    }
    
    method.vector <- rep(method, each = length(finalModel)/length(method))
    perf.metrics <- mapply(predicted,
                           FUN = function(x, y) perf.calc(x, lev = grp.levs, model = y),
                           y = method.vector, SIMPLIFY = FALSE)                 
    
    ## for classification, add the cell counts
    #require(plyr)
    
    #cells <- lapply(predicted,
    #                function(x) conf.matrix(x$pred, x$obs))
    cells <- lapply(predicted,
                    function(x) flatTable(x$pred, x$obs))
    for(ind in seq(along = cells)){
      perf.metrics[[ind]] <- c(perf.metrics[[ind]], cells[[ind]])
    } 
    
    final.metrics <- do.call("rbind", perf.metrics)
    
}

