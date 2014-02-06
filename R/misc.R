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

#' @title Feature Consistency Table
#' @description Extracts and sorts the features identified for a given method.
#' @param features A fs.stability fitted object
#' @param method Algorithm of interest 
#' Available options are \code{"plsda"} (Partial Least Squares Discriminant Analysis),
#'  \code{"rf"} (Random Forest), \code{"gbm"} (Gradient Boosting Machine),
#'  \code{"svm"} (Support Vector Machines), \code{"glmnet"} (Elastic-net Generalized Linear Model),
#'  and \code{"pam"} (Prediction Analysis of Microarrays)
#' @return A data frame containing:
#'  \item{features}{Features identified by model}
#'  \item{consistency}{Number of iterations feature was identified}
#'  \item{frequency}{Frequency of iterations the feature was identified}
#' @author Charles Determan Jr
#' @export


feature.table <- function(features, method = NULL){
  if(is.null(method)){
    stop("\n Error: must specify which algorithm to extract")}else{
      feats <- as.data.frame(lapply(features$features[method],'[[', "features"))
      k <- ncol(feats)
      table.dat <- sort(table(unlist(feats)), decreasing = T)
      out <- as.matrix(table.dat)
      out <- cbind(as.character(rownames(out)), as.numeric(out[,1]), round(as.numeric(out[,1])/k, 3))
      colnames(out) <- c("features", "consistency", "frequency")
      out <- data.frame(out)
      out
    }
}

# caret:::createFolds function

# Confusion Matrix generation
# import caret:::flatTable
#conf.matrix <- function(pred, obs)

# caret:::byComplexity
byComplexity2 <- function(x, model)
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
# caret:::best

# May not include this option
# caret:::oneSE

# caret:::tolerance


