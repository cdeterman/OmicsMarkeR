
#' @title Model Parameters and Properties
#' @description Provides a list of the models with their respective parameters 
#' and properties.
#' @param method A vector of strings listing the models to be returned
#' @return Returns a dataframe of the following components:
#' @return method A vector of strings listing models returned
#' @return parameter A vector of possible parameters to be optimized
#' @return label A vector of the names for each possible parameter
#' @return seq A logical indicator if the parameter is sequential in the 
#' model (i.e. if model is able to fit all 'lower' parameters simultaneously)
#' @example inst/examples/params.R
#' @export

params <- function(method = NULL)
{
    methods <- c(
        ## gbm
        'gbm', 'gbm', 'gbm',
        ## glmnet
        'glmnet', 'glmnet', 
        ## pam
        'pam', 
        ## plsda
        'plsda',
        ## rf
        'rf',
        ## svmLinear
        'svm')
    
    parameters <- c(
        ## gbm
        'n.trees', 'interaction.depth', 'shrinkage',
        ## glmnet
        'lambda', 'alpha', 
        ## pam
        'threshold', 
        ## plsda
        'ncomp',
        ## rf
        'mtry',
        ## svmLinear
        'C')
    
    
    labels <- c(
        ## gbm
        '#Trees',
        'Interaction Depth',
        'Learning Rate',
        ## glmnet
        'Regularization Parameter',
        'Mixing Percentage',
        ## pam
        'Shrinkage Threshold',
        ## plsda
        '#Components',
        ## rf
        '#Randomly Selected Predictors',
        ## svmLinear
        'C'
    )
    
    allone <- c(
        ## gbm
        TRUE, FALSE, FALSE,
        ## glmnet
        TRUE, FALSE, 
        ## pam
        TRUE, 
        ## pls
        TRUE,
        ## rf
        FALSE,
        ## svmLinear
        FALSE
    )
    
    params <- data.frame(method = methods,
                         parameter = parameters,
                         label = labels,
                         seq = allone)
    
    method <- tolower(method)
    
    if(!is.null(method))
    {
        if(!any(method %in% params$method)) stop("value of algorithm unknown")
        
        for (i in 1:length(method)){
            tmp <- list(params[params$method %in% method[i],])
            if(i == 1){
                out <- tmp
            }else{
                out <- c(out, tmp)
            }
        }
        names(out) <- method
        
    } else out <- params  
    
    out
}
