

verify <- function(x, y, method, na.rm = na.rm){  
    
    # check if x is matrix
    assert_is_matrix(x)
    # check if x is numeric
    assert_is_numeric(x)
    
    # check if y is factor
    assert_is_factor(y)
    # no missing values in y
    assert_all_are_not_na(y)
    
    # method must be character
    assert_is_character(method)
    
    # set method to lowercase to avoid user frustration
    method <- tolower(method)
    
    if(!all(method %in% c("plsda","glmnet","gbm", "pam", "rf", "svm"))){
        stop("Selected method is not plsda, glmnet, gbm, pam, 
                    rf or svm!")
    }
    
    # check lengths of x and y
    if (nrow(x) != length(y))
        stop("\n Error: 'variables' and 'group' have different lengths")
    # no missing values allowed when na.rm=FALSE
    if (!na.rm) {
        if (length(complete.cases(x)) != nrow(x))
            stop("\n Error: no missing values allowed in 'variables'")    
    }
}


verify_fs <- function(f, stability.metric,
                      model.features, no.fs)
{
    assert_is_character(stability.metric)
    assert_is_logical(model.features)
    assert_is_logical(no.fs)
    
    # added no.fs (i.e. no feature selection to reuse in fit.only.model 
    # function)
    if(!no.fs){
        # kuncheva is a special case stability metric
        if(stability.metric == "kuncheva"){
            if(is.null(f)){
                stop("Kuncheva requires same number of features 
                     in subsets. \nYou must specifiy number of features to 
                     be selected")    
            }
            if(model.features){
                stop("Kuncheva requires same number of features 
                     in subsets. \nYou must specifiy number of features to 
                     be selected (f).")
            }
            }
        
        # full rank corrleation doesn't allow for feature subsets
        if(model.features & !is.null(f)){
            stop("If running model defined number of features you can not 
                 specifiy the number of features.  \n")
        }
        # f cannot be set if doing rank corrleation metric
        if(stability.metric %in% c("spearman", "canberra") & !is.null(f)){
            warning("Rank Correlation doesn't allow for feature subsets.\n
                    'f' has been set to NULL")
            f <- NULL
        } 
        # make sure 'f' is set when using a feature subset metric
        if(stability.metric %in% c("jaccard","sorenson","ochiai",
                                   "pof","kuncheva") 
           & is.null(f) & !model.features){
            stop(paste("Stability metric", " '", stability.metric, "' ", 
                       "requires 'f' to be set.\n", sep = ""))
        }  
    }
    return(f)
}