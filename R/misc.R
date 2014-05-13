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


# slight modification of the 'caret' function 'byComplexity' to suit this program
byComplexity2 <- function(x, model)
{
  # must be dataframe to access components easily with '$'
  if(!is.data.frame(x)) x <- as.data.frame(x)
  
  switch(tolower(model),
         gbm =
{
  # As Max Kuhn (caret creator) has stated:
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



# These functions or duplicates from the 'caret' package
# They have been copied as they are low-level functions that are not exported
# However, they prove useful for this program as it has a very similar structure
MeanSD <- function(x, exclude = NULL)
{
  if(!is.null(exclude)) x <- x[, !(colnames(x) %in% exclude), drop = FALSE]
  out <- c(colMeans(x, na.rm = TRUE), sapply(x, sd, na.rm = TRUE))
  names(out)[-(1:ncol(x))] <- paste(names(out)[-(1:ncol(x))], "SD", sep = "")
  out
}

expandParameters <- function(fixed, seq)
{
  if(is.null(seq)) return(fixed)
  
  isSeq <- names(fixed) %in% names(seq)
  out <- fixed
  for(i in 1:nrow(seq))
  {
    tmp <- fixed
    tmp[,isSeq] <- seq[i,]
    out <- rbind(out, tmp)
  }
  out
}

flatTable <- function(pred, obs)
{
  cells <- as.vector(table(pred, obs))
  if(length(cells) == 0) cells <- rep(NA, length(levels(obs))^2)
  names(cells) <- paste(".cell", seq(along= cells), sep = "")
  cells
}

# Options for selecting optimal model by caret
# currently incorporating caret:::best
# Currently don't support oneSE or tolerance
# Dependent upon user needs
