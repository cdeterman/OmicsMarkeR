# sequester appropriate parameters

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


