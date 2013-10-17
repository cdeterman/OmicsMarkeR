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
                       full = allone)

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
