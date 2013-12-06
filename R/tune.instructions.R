
#' @title Model Optimization Instructions
#' @description Provides directions for which parameters to loop over during tuning.  This becomes
#' important when certain models can access 'lower' parameters without running them independently.
#' @param method Vector of strings indicating which models will be fit
#' @param grid A list of parameters grids to be applied to the models
#' @return \item{modelInfo}{List of the following components}.
#' @return \itemize{
#'  \item{scheme: String dictating which looping scheme to apply}
#'  \item{loop: Dataframe of parameters to loop through for each model}
#'  \item{model: Information regarding parameters of specific model}
#'  \item{constant: Names of the 'loop' dataframe components}
#'  \item{vary: Indication of parameters that vary and can access recursively}}
#' @author Charles E. Determan Jr.
#' @import DiscriMiner
#' @import randomForest
#' @import caret
#' @import e1071
#' @import gbm
#' @import pamr
#' @import glmnet
# ' @export

tune.instructions <- function(method, grid)
{
  modelInfo <- params(method)
  
  ## In a very similar fashion to the caret package, some models have parameters where 
  ## several different models can be derived from one R object. For example, in plsDA models 
  ## you can fit a model with 15 components and get predictions for any mode with <= 15 components 
  ## from the same object.
  
  ## if we don't have any of these types of parameters we can use a typical looping strategy
  ## i.e. scheme = "basic"
  
  for(i in 1:length(modelInfo)){
    if(!any(modelInfo[[i]]$seq)){
      modelInfo[[i]] <- 
        list( 
          scheme = "basic",
          loop = grid[[i]], 
          seqParam = NULL, 
          model = modelInfo[[i]], 
          constant = names(grid[[i]]), 
          vary = NULL)
    }
  }
  
  for(i in 1:length(grid)){
    if(any(modelInfo[[i]]$seq)){
      paramVary <- unlist(lapply(grid[[i]], function(u) length(unique(u)) > 1))
      paramVary <- data.frame(
        parameter = substring(names(paramVary), 2),
        column = names(paramVary),
        varies = paramVary)
      
      modelInfo[[i]] <- merge(modelInfo[[i]], paramVary)
      modelInfo[[i]]$varyingSeq <- modelInfo[[i]]$varies & modelInfo[[i]]$seq
      
      scheme <- if(any(modelInfo[[i]]$varyingSeq)) "seq" else "basic" 
      
      if(scheme == "seq")
      {
        constant <- as.character(modelInfo[[i]]$column)[!modelInfo[[i]]$varyingSeq]
        vary <- as.character(modelInfo[[i]]$column)[modelInfo[[i]]$varyingSeq] 
        
        ## The data frame loop is the combination(s) of tuning parameters that we will
        ## be looping over. For each combination in loop, the list seqParam will provide the
        ## value(s) of the sequential parameter that should be evaluated for the same R model
        ## object      
        
        switch(method[[i]],
               plsda = 
                 {
                   grid[[i]] <- grid[[i]][order(grid[[i]]$.ncomp, decreasing = TRUE),, drop = FALSE]
                   loop <- grid[[i]][1,,drop = FALSE]
                   seqParam <- list(grid[[i]][-1,,drop = FALSE])
                 },
               gbm =
                 { 
                   loop <- aggregate(
                     grid[[i]]$.n.trees, 
                     list(
                       .interaction.depth = grid[[i]]$.interaction.depth, 
                       .shrinkage = grid[[i]]$.shrinkage), max)
                   
                   names(loop)[3] <- ".n.trees"
                   seqParam <- vector(mode = "list", length = nrow(loop))
                   for(k in seq(along = loop$.n.trees))
                   {
                     index <- which(
                       grid[[i]]$.interaction.depth == loop$.interaction.depth[k] & 
                         grid[[i]]$.shrinkage == loop$.shrinkage[k])
                     subTrees <- grid[[i]][index, ".n.trees"] 
                     seqParam[[k]] <- data.frame(.n.trees = subTrees[subTrees != loop$.n.trees[k]])
                   }         
                 },
               pam = 
                 {
                   grid[[i]] <- grid[[i]][order(grid[[i]]$.threshold, decreasing = TRUE),, drop = FALSE]
                   loop <- grid[[i]][1,,drop = FALSE]
                   seqParam <- list(grid[[i]][-1,,drop = FALSE])
                 },
               glmnet =
                 {  
                   uniqueAlpha <- unique(grid[[i]]$.alpha)
                   
                   loop <- data.frame(.alpha = uniqueAlpha)
                   loop$.lambda <- NA
                   
                   seqParam <- vector(mode = "list", length = length(uniqueAlpha))
                   
                   for(k in seq(along = uniqueAlpha))
                   {
                     seqParam[[k]] <- data.frame(.lambda = subset(grid[[i]], subset = loop$.alpha == uniqueAlpha[k])$.lambda)
                   } 
                 }
        )
        
        modelInfo[[i]] <- list(scheme = "seq", 
                               loop = loop, 
                               seqParam = seqParam, 
                               model = modelInfo[[i]], 
                               constant = constant, 
                               vary = vary)
      } else {
        modelInfo[[i]] <- list(scheme = "basic", 
                               loop = grid[[i]], 
                               seqParam = NULL, 
                               model = modelInfo[[i]], 
                               constant = names(grid[[i]]), 
                               vary = NULL)
      }
    }else{
      next
    }
  }
  modelInfo
}