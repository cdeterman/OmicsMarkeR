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









