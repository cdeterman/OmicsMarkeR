

verify <- function(x, y, method, f, stability.metric, model.features, na.rm = na.rm){  
  
  # set method to lowercase to avoid user frustration
  method <- tolower(method)
  
  if(!all(method %in% c("plsda","glmnet","gbm", "pam", "rf", "svm")))
    simpleError("Selected method is not plsda, glmnet, gbm, pam, rf or svm!")
  
  # kuncheva is a special case stability metric
  if(stability.metric == "kuncheva"){
    if(is.null(f)){
      stop("\n Error: Kuncheva requires same number of features in subsets.\nYou must specifiy number of features to be selected")    
    }
    if(!model.features){
      stop("\n Error: Kuncheva requires same number of features in subsets. \nYou must specifiy number of features to be selected")
    }
  }
  
  # full rank corrleation doesn't allow for feature subsets
  if(model.features & !is.null(f)){
    stop("If running model defined number of features you can not specifiy the number of features.  \n")
  }
  # f cannot be set if doing rank corrleation metric
  if(stability.metric %in% c("spearman", "canberra")){
    warning("Rank Correlation doesn't allow for feature subsets.\n'f' has been set to NULL")
    f <- NULL
  } 
  # make sure 'f' is set when using a feature subset metric
  if(stability.metric %in% c("jaccard","sorenson","ochiai","pof","kuncheva") & is.null(f) & !model.features){
    stop(paste("Stability metric", " '", stability.metric, "' ", "requires 'f' to be set.\n", sep = ""))
  }  
  # x matrix or data.frame
  if (is.null(dim(x))) {
    stop("\n Error: 'variables' is not a matrix or dataframe")
  }
  # check lengths of x and y
  if (nrow(x) != length(y))
    stop("\n Error: 'variables' and 'group' have different lengths")
  # no missing values allowed when na.rm=FALSE
  if (!na.rm) {
    if (length(complete.cases(x)) != nrow(x))
      stop("\n Error: no missing values allowed in 'variables'")    
  }
  # y vector or factor
  if (!is.vector(y) && !is.factor(y))
    stop("\n Error: 'group' must be a factor")
  # make sure y is a factor
  if (!is.factor(y)) y = as.factor(y)
  # no missing values in y
  if (any(!is.finite(y)))
    stop("\n Error: no missing values allowed in 'group'")
  # quantitative data
  # make sure is matrix
  if (!is.matrix(x)) x <- as.matrix(x)
  # only numeric values
  if (!is.numeric(x))
    stop("\n Error: 'variables' must contain only numeric values")
  
  # verified inputs
  if (is.null(colnames(x))) colnames(x) = paste(rep("X",ncol(x)),seq_len(ncol(x)), sep='') 
  if (is.null(rownames(x))) rownames(x) = 1:nrow(x)
  
  list(X = x, Y = y, method = method, f = f)
}
