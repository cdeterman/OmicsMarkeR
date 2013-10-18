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

# caret:::best

# May not include this option
# caret:::oneSE

# caret:::tolerance









