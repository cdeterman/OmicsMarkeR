
#' @title Model List
#' @description Provide a list of currently implemented methods 
#' for OmicsMarkeR.
#' @return A data.frame containing:
#' \item{methods}{The abbreviated code for the method}
#' \item{description}{Full name of the method}
#' @author Charles Determan Jr.
#' @examples modelList()
#' @export
modelList <- function(){
    methods <- c(
        "plsda",
        "svm",
        "rf",
        "gbm",
        "glmnet",
        "pam")
    
    description <- c(
        "Partial Least Squares Discriminant Analysis",
        "Support Vector Machines",
        "Random Forest",
        "Gradient Boosting Machines",
        "Elastic-net regularized models",
        "Predication Analysis of Microarrays")
    
    out <- data.frame(methods = methods,
                      description = description)
    return(out)
}