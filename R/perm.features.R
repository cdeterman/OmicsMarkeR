### Monte Carlo Permutation Testing for Feature Significance

## Function to extract arguments from previously fit model
#' @title Argument extractor
#' @description Extract arguments from previously fs.stability models
#' @param fs.model Previously fit fs.stability model
#' @param method Which model to extract from
#' @return args List of model arguments

extract.args <- function(fs.model, method)
{
  ind <- which(names(fs.model$performance) == method)
  metrics <- data.frame(fs.model$performance[[ind]])
  
  args <- switch(method,
                 plsda = {out <- list(ncomp = metrics$ncomp)},
                 gbm = {out <- c(metrics[,colnames(metrics) %in% c("interaction.depth", "n.trees", "shrinkage")])},
                 rf = {out <- list(mtry = metrics$mtry)},
                 svm = {out <- list(C = metrics$C)},
                 pam = {out <- list(threshold = metrics$threshold)},
                 glmnet = {out <- c(metrics[colnames(metrics) %in% c("alpha", "lambda")])}
  )
  return(args)
}



#' @title Feature Selection via Monte Carlo Permutation
#' @description Applies Monte Carlo permutations to user specified models.  The user can either use the results
#' from \code{fs.stability} or provide specified model parameters.
#' @param fs.model Object containing results from \code{fs.stability}
#' @param X A scaled matrix or dataframe containing numeric values of each feature
#' @param Y A factor vector containing group membership of samples
#' @param method A vector listing models to be fit.
#' Available options are \code{"plsda"} (Partial Least Squares Discriminant Analysis),
#'  \code{"rf"} (Random Forest), \code{"gbm"} (Gradient Boosting Machine),
#'  \code{"svm"} (Support Vector Machines), \code{"glmnet"} (Elastic-net Generalized Linear Model),
#'  and \code{"pam"} (Prediction Analysis of Microarrays)
#' @param sig.level Desired significance level for features, default \code{sig.level = .05}
#' @param nperm Number of permutations, default \code{nperm = 10}
#' @param allowParallel Logical argument dictating if parallel processing is allowed via foreach package.
#' Default \code{allowParallel = FALSE}
#' @return \item{sig.level}{User-specified significance level}
#' @return \item{num.sig.features}{Number of significant features}
#' @return \item{sig.features}{Dataframe of significant features}
#' @author Charles Determan Jr.
#' @references Wongravee K., et. al. (2009) \emph{Monte-Carlo methods for determining optimal number
#' of significant variables.  Application to mouse urinary profiles}. Metabolomics 5:387-406.
#' @import DiscriMiner
#' @import randomForest
#' @import e1071
#' @import gbm
#' @import pamr
#' @import glmnet
#' @importFrom permute shuffle
#' @export

# user-defined -> sig.level, model, nperm
perm.features <- function(fs.model = NULL, X, Y, method, sig.level = .05, nperm = 10, allowParallel = FALSE, ...)
  {
  `%op%` <- if(allowParallel){
    `%dopar%` 
  }else{
    `%do%`
  } 
  
  theDots <- list(...)
  if(is.null(fs.model) & length(theDots) == 0){
    stop("Error: you must either provide fitted model from fs.stability or the parameters for the desired model")
  }
  
  obsLevels <- levels(as.factor(Y))
  
  data <- as.data.frame(cbind(X, .classes = Y))
  
  ## pam and will crash if there is a resample with <2 observations
  ## in a class. We will detect this and remove those classes.
  if(method == "pam")
  {
    yDist <- table(data$.classes)
    if(any(yDist < 2))
    {
      smallClasses <- names(yDist[yDist < 2])
      data <- data[!(data$.classes %in% smallClasses),]
    }
  }
  
  ## Factor the class labels
  data$.classes <- factor(as.character(data$.classes), levels = obsLevels)
  xNames <- names(data)[!(names(data) %in% ".classes")]
  
  trainX <- data[,!(names(data) %in% ".classes"), drop = FALSE]
  trainY <- data[,".classes"]
  
  if(method == "gbm" & length(obsLevels) == 2)  numClasses <- ifelse(data$.classes == obsLevels[1], 1, 0)
  
  if(!is.null(fs.model)){
    args <- extract.args(fs.model, method)
  }else{
    if(!is.null(theDots) & length(theDots) != 0){
      args <- theDots
    }
  }
  
  if(!is.null(theDots) & length(theDots) != 0){
    if(names(theDots) %in% names(args)){
      arg.ind <- which(names(theDots) %in% names(args))
      args <- theDots[arg.ind]
      theDots <- theDots[-arg.ind]
    }  
  }
  
  args <- data.frame(do.call("cbind", args))
  
  # extract new values of variables
  N <- nrow(trainX)
  # create permutations
  perms <- replicate(nperm, shuffle(N), simplify = F)
  perms <- append(list(seq(N)), perms)
  
  scores <- foreach(p = seq.int(nperm+1),
                    .packages = c("OmicsMarkeR", "permute", "foreach", "randomForest", "e1071", "DiscriMiner","gbm","pamr","glmnet"),
                    .verbose = FALSE,
                    .errorhandling = "stop") %op% 
    {
    
    features <- switch(method,
                       plsda ={
                         mod <- plsDA(trainX, 
                                      trainY,
                                      autosel=F,
                                      validation = NULL,
                                      comps = if(args$ncomp > 1) args$ncomp else 2,
                                      cv ="none")
                         
                         if(args < 2){
                           l <- 1
                         }else{
                           l <- ncol(mod$VIP)  
                         }
                         out <- mod$VIP[,l]
                       },
                       
                       gbm ={
                         gbm.args <- c("w", "var.monotone", "n.minobsinnode", 
                                       "bag.fraction", "var.names", "response.name", "group") 
                         theDots <- theDots[names(theDots) %in% gbm.args]
                         
                         if(ncol(trainX) < 50 | nrow(trainX) < 50){
                           if(is.null(theDots) | length(theDots) == 0){
                             if(nrow(trainX) < 30){
                               theDots <- list(n.minobsinnode = 2)
                             }else{
                               theDots <- list(n.minobsinnode = 5)  
                             }
                           }
                         }
                         # determine if binary or multiclass
                         gbmdist <- if(length(unique(trainY)) == 2){
                           "bernoulli"}else{
                             "multinomial"
                           }         
                         
                         # check gbm setup file to see if this is necessary
                         modY <- if(gbmdist != "multinomial") numClasses else trainY
                         
                         if(gbmdist != "multinomial"){
                           modY <- numClasses
                         }else{
                           modY <- trainY
                         }
                         
                         modArgs <- list(x = X,
                                         y = modY,
                                         interaction.depth = args$interaction.depth,
                                         n.trees = args$n.trees,
                                         shrinkage = args$shrinkage, 
                                         distribution = gbmdist,
                                         verbose = FALSE)
                         
                         if(length(theDots) > 0) modArgs <- c(modArgs, theDots)
                         
                         mod <- do.call("gbm.fit", modArgs)
                         
                         out <- relative.influence(mod, n.trees = args$n.trees)
                       },
                       
                       rf ={
                         rf.args <- c("maxnodes", "keep.forest", "keep.inbag")
                         theDots <- theDots[names(theDots) %in% rf.args]
                         
                         modArgs <- list(x = trainX,
                                         y = trainY,
                                         importance = TRUE,
                                         mtry = args$mtry,
                                         ntree=round.multiple(sqrt(ncol(trainX)), target = 50)
                         )
                         
                         if(length(theDots) > 0) modArgs <- c(modArgs, theDots)
                         
                         mod <- do.call("randomForest", modArgs)
                         
                         # Mean Decrease in Accuracy metric (type=1)
                         # Gini Index (type=2)
                         out <- importance(mod, type = 1)
                       },                  
                       
                       svm ={
                         best.C <- c(args$C)
                         
                         if(nlevels(trainY) == 2){
                           out <- svmrfeFeatureRanking(trainX, trainY, best.C)
                         }else{
                           out <- svmrfeFeatureRankingForMulticlass(trainX, trainY, best.C)
                         }
                       },
                       
                       pam ={
                         pamr.args <- c("n.threshold", "threshold.scale", "scale.sd", "se.scale")
                         theDots <- theDots[names(theDots) %in% pamr.args]
                         #p=100
                         modArgs <- list(data = list(x = t(trainX), y = trainY[perms[[p]]], geneid = as.character(colnames(trainX))),
                                         threshold = args$threshold
                         )
                         mydata <- list(x = t(trainX), y = trainY[perms[[p]]], geneid = as.character(colnames(trainX)))
                         
                         
                         if(length(theDots) > 0) modArgs <- c(modArgs, theDots)
                         
                         mod <- do.call("pamr.train", modArgs)
                         if(mod$nonzero == 0){
                           zero <- replicate(length(obsLevels), rep(0, length(xNames)))
                           pam.features <- data.frame(factor(xNames), zero)
                           colnames(pam.features) <- c("id", paste(seq(length(obsLevels)), "score", sep = "-"))
                         }else{
                           pam.features <- data.frame(pamr.listgenes(mod,
                                                                     mydata,
                                                                     threshold = args$threshold)
                           )
                         }
                         
                        
                        pam.scores <- sapply(pam.features[,2:ncol(pam.features)], FUN = function(x) as.numeric(as.character(x)))
                         # check to make sure in matrix format
                         if(!is.matrix(pam.scores)){
                           pam.scores <- t(as.data.frame(pam.scores))
                         }
                         out <- rowSums(abs(pam.scores))
                         names(out) <- as.character(pam.features[,1])
                         out
                       },         
                       
                       glmnet ={
                         numLev <- if(is.character(trainY) | is.factor(trainY)) length(levels(trainY)) else NA
                         
                         glmnet.args <- c("offset", "nlambda", "weights", "standardize","intecept", 
                                          "dfmax", "pmax","exclude","penalty.factor","lower.limits",
                                          "upper.limits","maxit","standardize.response","type.multinomial")
                         theDots <- theDots[names(theDots) %in% glmnet.args]
                         
                         if(!is.null(theDots)){
                           if(all(names(theDots) != "family"))
                           {
                             if(!is.na(numLev))
                             {
                               fam <- ifelse(numLev > 2, "multinomial", "binomial")
                             } else stop("Error: levels of classes couldn't be determined for glmnet")
                             
                             if(is.null(theDots)){
                               theDots <- list(family = fam)
                             }else{
                               theDots$family <- fam
                             }
                           }
                         }else{
                           if(!is.na(numLev))
                           {
                             fam <- ifelse(numLev > 2, "multinomial", "binomial")
                           } else stop("Error: levels of classes couldn't be determined for glmnet")
                           
                           if(is.null(theDots)){
                             theDots <- list(family = fam)
                           }
                         }
                         
                         modelArgs <- c(
                           list(x = as.matrix(trainX), y = trainY, alpha = args$alpha), theDots)
                         
                         # fit the model
                         mod <- do.call("glmnet", modelArgs)
                         
                         # extract coefficients
                         full.coefs <- coef(mod, s = 0)
                         if(nlevels(trainY) > 2){
                           full.coefs <- lapply(full.coefs, FUN = function(x) abs(as.matrix(x[2:(ncol(trainX)+1),, drop = FALSE])))
                           coefs <- lapply(full.coefs, FUN = function(x) x[x[,1] >= 0,, drop = FALSE])
                           
                           coef.names <- lapply(coefs, row.names)
                           coefs <- unlist(coefs, use.names = TRUE)
                           names(coefs) <- unlist(coef.names)
                           
                           #length(unique(names(coefs)))
                           var.names <- unique(names(coefs))
                           
                           for(n in seq(along = unique(names(coefs)))){
                             ind <- which(names(coefs) == var.names[n])
                             uni <- sum(abs(coefs[ind]))
                             names(uni) <- var.names[n]
                             coefs <- coefs[-ind]
                             coefs <- c(coefs, uni)
                           }   
                           out <- coefs
                         }else{
                           out <- abs(as.matrix(full.coefs[2:(ncol(trainX)+1),, drop = FALSE]))
                         }
                         scs <- out[,1]
                         names(scs) <- rownames(out)
                         scs
                       }
    )
    features
  }
  
  # extract p-value (one-tailed)
  var.scores <- lapply(scores, "[")  
  miss <- lapply(var.scores, FUN = function(x) which(!xNames %in% names(x)))
  for(i in 1:(nperm+1)){
    if(length(miss[[i]]) >= 1){
      zero <- rep(0, length(miss[[i]]))
      names(zero) <- xNames[miss[[i]]]
      var.scores[[i]] <- c(var.scores[[i]], zero)
    } 
  }

  var.scores <- sapply(var.scores, FUN = function(x) x[sort(names(x))])
  perm.p.val=apply(var.scores, 1, FUN = function(x) sprintf("%.3f", round(sum(x[2:(nperm+1)] >= x[1])/nperm, digits = 3)))
  
  # return number of significant features and features themselves
  num.features <- sum(perm.p.val <= sig.level)
  features <- names(which(perm.p.val <= sig.level))
  p.vals <- perm.p.val[which(perm.p.val <= sig.level)]
  
  cat("\nFeature Permutation Results\n")
  cat(rep("-",30), sep="")
  cat("\n\n")
  cat(paste(num.features, "features were significant at significance level", sig.level, sep = " "))
  cat("\nIdentified features:\n")
  print(data.frame(features = features, p.val = p.vals))
  sig.features <- data.frame(features = features, p.val = p.vals)
  all.features <- data.frame(features = names(perm.p.val), p.val = c(perm.p.val))
  
  perm.results = list(sig.level = sig.level,
                      num.sig.features = num.features,
                      sig.features = sig.features,
                      all.features = all.features
                      )
  
  return(perm.results)
}




