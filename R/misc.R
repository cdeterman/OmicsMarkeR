########################################
### Miscellaneous Functions Required ###
########################################

# Round to nearest multiple of target number
round.multiple <- function(x, target, f = round) {
  out <- f(x / target) * target
  if(out == 0){
    out <- 50
  }
  out
}

# caret:::createFolds function
createFolds <- 
  function(y, k = 10, list = TRUE, returnTrain = FALSE)
  {
    
    if(is.numeric(y))
    {
      ## Group the numeric data based on their magnitudes
      ## and sample within those groups.
      
      ## When the number of samples is low, we may have
      ## issues further slicing the numeric data into
      ## groups. The number of groups will depend on the
      ## ratio of the number of folds to the sample size.
      ## At most, we will use quantiles. If the sample
      ## is too small, we just do regular unstratified
      ## CV
      cuts <- floor(length(y)/k)
      if(cuts < 2) cuts <- 2
      if(cuts > 5) cuts <- 5
      y <- cut(
        y, 
        unique(
          quantile(y,
                   probs =
                     seq(0, 1, length = cuts))), 
        include.lowest = TRUE)
    }
    
    
    if(k < length(y))
    {
      ## reset levels so that the possible levels and 
      ## the levels in the vector are the same
      y <- factor(as.character(y))
      numInClass <- table(y)
      foldVector <- vector(mode = "integer", length(y))
      
      ## For each class, balance the fold allocation as far 
      ## as possible, then resample the remainder.
      ## The final assignment of folds is also randomized. 
      for(i in 1:length(numInClass))
      {
        ## create a vector of integers from 1:k as many times as possible without 
        ## going over the number of samples in the class. Note that if the number 
        ## of samples in a class is less than k, nothing is producd here.
        seqVector <- rep(1:k, numInClass[i] %/% k)
        ## add enough random integers to get  length(seqVector) == numInClass[i]
        if(numInClass[i] %% k > 0) seqVector <- c(seqVector, sample(1:k, numInClass[i] %% k))
        ## shuffle the integers for fold assignment and assign to this classes's data
        foldVector[which(y == dimnames(numInClass)$y[i])] <- sample(seqVector)
      }
    } else foldVector <- seq(along = y)
    
    if(list)
    {
      out <- split(seq(along = y), foldVector)
      names(out) <- paste("Fold", gsub(" ", "0", format(seq(along = out))), sep = "")
      if(returnTrain) out <- lapply(out, function(data, y) y[-data], y = seq(along = y))
    } else out <- foldVector
    out
  }

# Confusion Matrix generation
# import caret:::flatTable
#conf.matrix <- function(pred, obs)

flatTable <- function(pred, obs)
{
  cells <- as.vector(table(pred, obs))
  if(length(cells) == 0) cells <- rep(NA, length(levels(obs))^2)
  names(cells) <- paste(".cell", seq(along= cells), sep = "")
  cells
}


# This function sorts the tuning parameter matrix from
# least complex models to most complex models
# caret:::byComplexity
byComplexity <- function(x, model)
  {
  switch(tolower(model),
         gbm =
           {
             # This is a toss-up, but the # trees probably adds
             # complexity faster than number of splits
             x[order(x$n.trees, x$interaction.depth, x$shrinkage),] 
           },
         rf =, plsda =, pam = 
           {
             x[order(x[,1]),]
           },
         svm =
           {
             x[order(x$C),]
           },
         glmnet = x[order(-x$lambda, x$alpha),],
         stop("no sorting routine for this model")
  )
}

## In these functions, x is the data fram of performance values and tuning parameters.

# caret:::best
best <- function(x,             # dataframe of performance values and model parameters
                 metric        # necessary?
                 #maximize       # logical value to say whether higher values equal better performance
                 )
{
  
  #bestIter <- if(maximize) which.max(x[,metric])
  #else which.min(x[,metric])  
  bestIter <- which.max(x[,metric])
  bestIter
}

# May not include this option
# caret:::oneSE
oneSE <- function(x, metric, num, maximize)
{
  index <- 1:nrow(x)
  
  if(!maximize)
  {
    bestIndex <- which.min(x[,metric])  
    perf <- x[bestIndex,metric] + (x[bestIndex,paste(metric, "SD", sep = "")])/sqrt(num)
    candidates <- index[x[, metric] <= perf]
    bestIter <- min(candidates)
  } else {
    bestIndex <- which.max(x[,metric])  
    perf <- x[bestIndex,metric] - (x[bestIndex,paste(metric, "SD", sep = "")])/sqrt(num)
    
    candidates <- index[x[, metric] >= perf]
    bestIter <- min(candidates)
  }
  bestIter
}

# caret:::tolerance
tolerance <- function(x,            # dataframe of performance values and model parameters
                      metric,       # necessary???
                      tol = 1.5    # tolerance, i.e. percent worse performance in exchange for simpler model
                      #maximize      # logical value to say whether higher values equal better performance
                      )
{
  index <- 1:nrow(x)
  
  #if(!maximize)
  #{
  #  best <- min(x[,metric])  
  #  perf <- (x[,metric] - best)/best * 100
  #  candidates <- index[perf < tol]
  #  bestIter <- min(candidates)
 # } else {
    best <- max(x[,metric])  
    perf <- (x[,metric] - best)/best * -100
    candidates <- index[perf < tol]
    bestIter <- min(candidates)
  #}
  bestIter
}


#############################
### Summarizing functions ###
#############################

# modified from caret:::postResample to include all model diagostics, applied only to classification problems
performance.stats <- function(pred, obs)
{
  isNA <- is.na(pred)
  pred <- pred[!isNA]
  obs <- obs[!isNA]
  
  if(length(obs) + length(pred) == 0)
  {
    out <- rep(NA, 2)
  } else {
    require(caTools)
    tmp.auc <- colMeans(colAUC(order(pred), obs, plotROC = FALSE, alg = "ROC"))
    pred <- factor(pred, levels = levels(obs))  
        
    tmp <- confusionMatrix(pred, obs)
    
    overall.stats <- c("Accuracy", "Kappa")
    byclass.stats <- c("Sensitivity", "Specificity", "Pos Pred Value", "Neg Pred Value")
    
    overall.index <- which(names(tmp$overall) %in% overall.stats)
    byclass.index <- which(names(tmp$byClass) %in% byclass.stats)
  
    out <- c(tmp$overall[overall.index], dim(tmp$table)[1], tmp$byClass[byclass.index], tmp.auc)
  }
  names(out) <- c("Accuracy", "Kappa", "Table", "Sensitivity", "Specificity", "Pos Pred Value", "Neg Pred Value", "ROC.AUC")
  out <- out[c("Accuracy", "Kappa", "ROC.AUC", "Sensitivity", "Specificity", "Pos Pred Value", "Neg Pred Value", "Table")]
  
  if(any(is.nan(out))) out[is.nan(out)] <- NA
  out
}


# modified from caret:::defaultSummary and caret:::twoClassSummary
perf.calc <- function(data,          
                      lev = NULL, 
                      model = NULL)     
{
  if(is.character(data$obs)){
    data$obs <- factor(data$obs, levels = lev)
    data$pred <- factor(data$pred, levels = lev)
  } 
  out <- performance.stats(data[,"pred"], data[,"obs"])
  out <- out[!names(out) == "Table"]
  out
}







