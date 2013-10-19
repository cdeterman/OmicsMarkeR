
#' @title Sequester Additional Parameters
#' @description When the user provides additional arguments to either \code{fs.stability} or 
#' \code{fs.ensembl.stability} this function will extract the parameters to be fit if optimization 
#' is not used i.e. \code{optimize = FALSE}.
#' @param theDots List of additional arguments
#' @param method Vector of strings listing models to be fit
#' @return Returns a list of the following elements
#' @return \item{parameters}{The parameters that will be fit to models}
#' @return \item{pnames}{The names of the specific parameters}

sequester <- function(theDots, method){

  Value <- vector("list", length(method))
  pnames <- vector()
  names(Value) <- method
  
  for(m in seq(along = method)){
    vals <- switch(tolower(method[m]),
                   plsda =
                     {
                       x <- which(names(theDots) == ".ncomp")
                       if(length(x) < 1){
                         stop("Error: If not autotuning or providing a grid, PLSDA requires you to specify 'ncomp'")
                       }
                       Value[[m]] <- theDots[x]
                       pnames <- c(pnames, "ncomp")
                     },
                   
                   gbm =
                     {
                       x <- which(names(theDots) %in% c(".n.trees", ".interaction.depth", ".shrinkage"))
                       if(length(x) < 3){
                         stop("Error: If not autotuning or providing a grid, GBM requires you to specify 'n.trees', 'interaction.depth', and 'shrinkage' values")
                       }
                       Value[[m]] <- theDots[x]
                       pnames <- c(pnames, c("n.trees", "interaction.depth", "shrinkage"))
                     },
                   
                   rf =
                     {
                       x <- which(names(theDots) == ".mtry")
                       if(length(x) < 1){
                         stop("Error: If not autotuning or providing a grid, RandomForest requires you to specify 'mtry'")
                       }
                       Value[[m]] <- theDots[x]
                       pnames <- c(pnames, "mtry")
                     },
                   
                   svm = 
                     {
                       x <- which(names(theDots) == ".C")
                       if(length(x) < 1){
                         stop("Error: If not autotuning or providing a grid, SVM requires you to specify 'C'")
                       }
                       Value[[m]] <- theDots[x]
                       pnames <- c(pnames, "C")
                     },
                   
                   glmnet = 
                     {
                       x <- which(names(theDots) %in% c(".lambda", ".alpha"))
                       if(length(x) < 1){
                         stop("Error: If not autotuning or providing a grid, glmnet requires you to specify 'lambda' and 'alpha'")
                       }
                       Value[[m]] <- theDots[x]
                       pnames <- c(pnames, c("lambda", "alpha"))
                     },
                   
                   pam = 
                     {
                       x <- which(names(theDots) == ".threshold")
                       if(length(x) < 1){
                         stop("Error: If not autotuning or providing a grid, PAM requires you to specify 'threshold'")
                       }
                       Value[[m]] <- theDots[x]
                       pnames <- c(pnames, "threshold")
                     }
    )
  }
 out <- list(parameters = Value,
             pnames = pnames)
  out
  
}


