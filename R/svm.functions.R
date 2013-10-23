
#' @title SVM Recursive Feature Extraction (Binary)
#' @description This conducts feature selection for Support Vector Machines models via recursive 
#' feature extraction.  This returns a vector of the features in x ordered by relevance.
#' The first item of the vector has the index of the feature which is more relevant to perform the
#' classification and the last item of the vector has the feature which is less relevant.
#' This function is specific to Binary classification problems,
#' @param x A matrix where each column represents a feature and each row represents a sample
#' @param y A vector of labels corresponding to each sample's group membership
#' @param c A numeric value corresponding to the 'cost' applied during the svm model fitting.
#' This can be selected by the user if using this function directly or is done internally.
#' @param perc.rem A numeric value indicating the percent of features removed during each iteration.
#' Default \code{perc.rem = 10}.
#' @return Vector of features ranked from most important to least important.
#' @references Guyon I. et. al. (2010) \emph{Gene Selection for Cancer Classification using Support Vector Machines}. 
#' Machine Learning 46 389-422. 
#' 
#' @seealso \code{\link{svmrfeFeatureRankingForMulticlass}}
#' @import e1071
#' @export

svmrfeFeatureRanking = function(x, y, c, perc.rem=10)
  {
  n = ncol(x)
  survivingFeaturesIndexes = seq(n)
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
    r <- ceiling(length(survivingFeaturesIndexes)/perc.rem)
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

#' @title SVM Multiclass Weights Ranking
#' @description This calculates feature weights for multiclass Support Vector Machine (SVM) problems
#' @param model A fitted SVM model of multiclass
#' @return Vector of feature weights
#' @references Guyon I. et. al. (2010) \emph{Gene Selection for Cancer Classification using Support Vector Machines}. 
#' Machine Learning 46 389-422. 

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

#' @title SVM Recursive Feature Extraction (Multiclass)
#' @description This conducts feature selection for Support Vector Machines models via recursive 
#' feature extraction.  This returns a vector of the features in x ordered by relevance.
#' The first item of the vector has the index of the feature which is more relevant to perform the
#' classification and the last item of the vector has the feature which is less relevant.
#' This function is specific to Binary classification problems,
#' @param x A matrix where each column represents a feature and each row represents a sample
#' @param y A vector of labels corresponding to each sample's group membership
#' @param c A numeric value corresponding to the 'cost' applied during the svm model fitting.
#' This can be selected by the user if using this function directly or is done internally.
#' @param perc.rem A numeric value indicating the percent of features removed during each iteration.
#' Default \code{perc.rem = 10}.
#' @return Vector of features ranked from most important to least important.
#' @references Guyon I. et. al. (2010) \emph{Gene Selection for Cancer Classification using Support Vector Machines}. 
#' Machine Learning 46 389-422. 
#' @seealso \code{\link{svmrfeFeatureRanking}}
#' @export

svmrfeFeatureRankingForMulticlass = function(x,y,c, perc.rem = 10){
  n = ncol(x)
  survivingFeaturesIndexes = seq(n)
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
    
    ## New to determine prec.rem of remaining features
    #round to make sure a defined number of features to remove (ceiling rounds up)
    r <- ceiling(length(survivingFeaturesIndexes)/perc.rem)
    s <- length(survivingFeaturesIndexes)
    
    #update feature ranked list
    #(featureRankedList[rankedFeatureIndex] = survivingFeaturesIndexes[ranking[1]])
    ## New update feature ranked list
    featureRankedList[rev((s-r+1):s)] <- survivingFeaturesIndexes[ranking[1:r]]
    
    ## Old to remove just 1
    #rankedFeatureIndex = rankedFeatureIndex - 1
    
    ## New to remove perc.rem
    rankedFeatureIndex <- rankedFeatureIndex - r
    
    #eliminate the feature with smallest ranking criterion
    survivingFeaturesIndexes = survivingFeaturesIndexes[-ranking[1:r]]
    cat(length(survivingFeaturesIndexes),"\n")
  }
}

