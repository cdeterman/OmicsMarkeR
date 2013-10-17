#################################################
### Support Vector Machines Feature Selection ###
#################################################

# working on facilitating remove 10%, currently just provides the last rank 300 times
#x <- subsample
#y <- subgroup
#rm(x,y)

## SVM Recursive feature extraction algorithms
## Based upon "Gene Selection for Cancer Classification using Support Vector Machines" by Isabelle Guyon 2002

##  The function has two inputs:
#   x : a matrix where each column represents a feature and each row represents a sample
#   y : a vector of the labels corresponding to each sample

#   The function returns a vector of the features in x ordered by relevance. The first item of
#   the vector has the index of the feature which is more relevant to perform the
#   classification and the last item of the vector has the feature which is less relevant.


svmrfeFeatureRanking = function(x, y, c)
  {
  n = ncol(x)
  survivingFeaturesIndexes = seq(1:n)
  featureRankedList = vector(length=n)
  rankedFeatureIndex = n
  while(length(survivingFeaturesIndexes)>0){
    #train the support vector machine
    svmModel <- svm(x[, survivingFeaturesIndexes], y, cost = c, cachesize=500,
                    scale=F, type="C-classification", kernel="linear" )
    #compute the weight vector
    w <- t(svmModel$coefs)%*%svmModel$SV
    #compute ranking criteria
    rankingCriteria <- w * w
    #rank the features
    ranking <- sort(rankingCriteria, index.return = TRUE)$ix

    ## New to determine 10% of remaining features
    #round to make sure a defined number of features to remove (ceiling rounds up)
    r <- ceiling(length(survivingFeaturesIndexes)/10)
    s <- length(survivingFeaturesIndexes)
    
    ## Old update feature ranked list
    #featureRankedList[rankedFeatureIndex] <- survivingFeaturesIndexes[ranking[1]]
    
    ## New update feature ranked list
    #lowestRankIndex <- rev((n-r-1):n)
    featureRankedList[rev((s-r+1):s)] <- survivingFeaturesIndexes[ranking[1:r]]
    
    ## Old to remove just 1
    #rankedFeatureIndex = rankedFeatureIndex - 1
    
    ## New to remove 10%
    rankedFeatureIndex <- rankedFeatureIndex - r
    
    #eliminate the feature with smallest ranking criterion
    survivingFeaturesIndexes <- survivingFeaturesIndexes[-ranking[1:r]]
  }
  return (featureRankedList)
}



## SVM Recursive Feature Extraction Algorithm (multiclass)
svm.weights<-function(model){
  w=0
  if(model$nclasses==2){
    w=t(model$coefs)%*%model$SV
  }else{ #when we deal with OVO svm classification
    ## compute start-index
    start <- c(1, cumsum(model$nSV)+1)
    start <- start[-length(start)]
    calcw <- function (i,j) {
      ## ranges for class i and j:
      ri <- start[i] : (start[i] + model$nSV[i] - 1)
      rj <- start[j] : (start[j] + model$nSV[j] - 1)
      ## coefs for (i,j):
      coef1 <- model$coefs[ri, j-1]
      coef2 <- model$coefs[rj, i]
      ## return w values:
      w=t(coef1)%*%model$SV[ri,]+t(coef2)%*%model$SV[rj,]
      return(w)
    }
    W=NULL
    for (i in 1 : (model$nclasses - 1)){
      for (j in (i + 1) : model$nclasses){
        wi=calcw(i,j)
        W=rbind(W,wi)
      }
    }
    w=W
  }
  return(w)
}

svmrfeFeatureRankingForMulticlass = function(x,y,c){
  n = ncol(x)
  survivingFeaturesIndexes = seq(1:n)
  featureRankedList = vector(length=n)
  rankedFeatureIndex = n
  while(length(survivingFeaturesIndexes)>0){
    #train the support vector machine
    svmModel = svm(x[, survivingFeaturesIndexes], y, cost = c, cachesize=500,
                   scale=F, type="C-classification", kernel="linear" )
    #compute the weight vector
    multiclassWeights = svm.weights(svmModel)
    #compute ranking criteria
    multiclassWeights = multiclassWeights * multiclassWeights
    rankingCriteria = 0
    for(i in 1:ncol(multiclassWeights))rankingCriteria[i] =
      mean(multiclassWeights[,i])
    #rank the features
    (ranking = sort(rankingCriteria, index.return = TRUE)$ix)
    #update feature ranked list
    (featureRankedList[rankedFeatureIndex] = survivingFeaturesIndexes[ranking[1]])
    rankedFeatureIndex = rankedFeatureIndex - 1
    #eliminate the feature with smallest ranking criterion
    (survivingFeaturesIndexes = survivingFeaturesIndexes[-ranking[1]])
    cat(length(survivingFeaturesIndexes),"\n")
  }
}
