### create grids

rfTune <- function(
  data,    # training data subset
  res
)
{
  require(randomForest)
  p <- dim(data)[2] - 1 
  c <- ceiling(p/50)
  
  # sequence of trees to try
  # both for high and 'very high' dimensional data
  if(p < 500 ){
    if(p < 100){
      treeSeq <- floor(seq(1, to = p, length = res))    
    }else{
      treeSeq <- floor(seq(10, to = p, length = c))
    }
  }else{
    treeSeq <- floor(2^seq(5, to = log(p, base = 2), length = c))
  } 
  
  # check if any of the numbers are repeated (i.e. repeating the same number of trees is inefficient)
  if(any(table(treeSeq) > 1))
  {
    treeSeq <- unique(treeSeq)
  }
  data.frame(.mtry = treeSeq)
}

## We fit an initial algorithm to get the range thresholds for these data and
## create the grid from these.
pamTune <- function(data,   # training data subset
                    res
)
{
  require(pamr)
  data <- data[complete.cases(data),,drop = FALSE]
  
  train.x <- data[!(names(data) %in% ".classes")]
  train.y <- data[,".classes"]
  
  
  initialThresh <- pamr.train(list(x=t(train.x), y=train.y))$threshold
  # remove 0.000 and highest thresholds
  initialThresh <- initialThresh[-c(1, length(initialThresh))]
  threshSeq <- data.frame(
    .threshold = seq(
      from = min(initialThresh),
      to = max(initialThresh), length = res))
  ## pamr.train prints out cv iterations without a line break
  #cat("\n")         
  threshSeq
}

# original caret function createGrid.R
denovo.grid <- function(data,       # training data of method being tuned
                        method,     # which algorithm to create grid
                        res         # resolution (i.e. how many/fine divisions of parameters)
                        
)
{
  for(i in 1:length(method)){
    algo <- tolower(method[i])
    tmp <- switch(algo,
                  # .ncomp doesn't need to be looped, can access any lesser components that max
                  plsda = expand.grid(
                    .ncomp = seq(1,res)),
                  
                  # alpha always tuned
                  # lambda tuning depends on user decision
                  #   if defined number of features (f), cannot be tuned
                  #   if doing full rank correlations (i.e. f = nc), cannot be tuned
                  #   if letting glmnet decide optimal features, then tune
                  glmnet = expand.grid(
                    .alpha = seq(0.1, 1, length = res),
                    .lambda = seq(.1, 3, length = 3 * res)),
                  
                  # n.trees doesn't need to be looped because can access any prior number from max
                  # modified .shrinkage from just '.1' to sequence of values
                  gbm = expand.grid(
                    .interaction.depth = seq(1, res),
                    .n.trees = floor((1:res) * 50),
                    .shrinkage = c(.1/seq(res))),
                  
                  rf = rfTune(data, res),
                  
                  pam = pamTune(data, res),
                  
                  svm = data.frame(.C = 2 ^((1:res) - 3)))
    
    if(i == 1){
      out <- list(tmp)
    }else{
      out <- append(out, list(tmp))
    }
  }
  names(out) <- tolower(method)
  out
}