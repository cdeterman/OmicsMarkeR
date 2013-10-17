
#####################################
### Stability/Similarity Measures ###
#####################################

# Spearman Rank Correlation Coefficient
spearman <- function(x,y){
  x <- as.vector(as.numeric(as.character(x)))
  y <- as.vector(as.numeric(as.character(y)))
  if(length(x) != length(y)){
    stop("\n Error: feature lengths must be same for Spearman Correlation")
  }
  N <- length(x)
  
  out <- 1 - 6*sum(((x - y)^2)/(N*((N^2)-1)))
  return(out)
}
  
## Canberra Distance
# formula from He 2010 (MR2) Jurman 2008
# higher value = greater distance apart (i.e. less similar)
canberra <- function(x,y){
  x <- as.vector(as.numeric(as.character(x)))
  y <- as.vector(as.numeric(as.character(y)))
  sum(abs(x-y)/(x+y))
}

## Jaccard Index
# formula in He 2010 (MS2) 
jaccard <- function(x,y){
  x <- as.vector(x)
  y <- as.vector(y)
  index <- length(intersect(x,y))/length(union(x, y))
  return(index)
}

## Dice-Sorensen's Index
# formula from He 2010 (MS3)
sorensen <- function(x,y){
  x <- as.vector(x)
  y <- as.vector(y)
  index <- 2*(length(intersect(x,y)))/(2*(length(intersect(x,y)))+length(setdiff(x,y))+length(setdiff(y,x)))
  return(index)
}

## Ochiai's Index
# formula from He 2010 (MS4)
ochiai <- function(x,y){
  x <- as.vector(x)
  y <- as.vector(y)
  index <- 2*(length(intersect(x,y)))/(sqrt(length(x)*length(y)))
  return(index)
}

## Percentage overlapping features
# formula from He 2010 (MS5)
pof <- function(x,y){
  x <- as.vector(x)
  y <- as.vector(y)
  index <- length(intersect(x,y))/length(x)
  return(index)
}

## Kuncheva Index
# formula from He 2010 (MS7)
# better description in Han 2012 - formula used below
# m/N = total number of features
# s/c = length of features - assumes equal length of x and y
# r = number of common elements in both signatures
# (r*N - s^2)/s*(N-s)
kuncheva <- function(x,             # list of features in first run
                     y,             # list of features in second run
                     num.features   # total number of features in entire dataset
                     )
  {
  x <- as.vector(x)
  y <- as.vector(y)
  if(length(x) != length(y)){
    stop("\n Error: Kuncheva requires same number of features in subset.\nYou must specifiy number of features to be selected")    
  }
  
  r <- length(intersect(x,y))
  d <- num.features
  k <- length(x)      # cardinality
  
  index <- (r - ((k^2)/d))/(k - ((k^2)/d))
  
  return(index)
}

#' @title Pairwise Stability Metrics
#' @description Conducts all pairwise comparisons of features selected following bootstrapping.
#' Also known as the data perturbation ensemble approach
#' @param features A matrix of selected features
#' @param stability.metric string indicating the type of stability metric.
#' Available options are \code{"jaccard"} (Jaccard Index/Tanimoto Distance),
#'  \code{"sorensen"} (Dice-Sorensen's Index), \code{"ochiai"} (Ochiai's Index),
#'  \code{"pof"} (Percent of Overlapping Features), \code{"kuncheva"} (Kuncheva's Stability Measures),
#'  \code{"spearman"} (Spearman Rank Correlation), and \code{"canberra"} (Canberra Distance)
#' @param orig.nc Original number of variables (Only applicable when using \code{"Kuncheva"})
#' @return A list is returned containing:
#' \item{comparisons}{Matrix of pairwise comparisons}
#' \item{overall}{The average of all pairwise comparisons}
#' @author Charles Determan Jr
#' @references He. Z. & Weichuan Y. (2010) \emph{Stable feature selection for biomarker discovery}. 
#' Computational Biology and Chemistry 34 215-225.
pairwise.stability <- 
  function(features, 
           stability.metric, 
           nc
           )
    {
    k <- ncol(features)
    comp <- matrix(0, nrow=k-1, ncol=k)
    for (i in 1:(k-1)){
      for (j in (i+1):k){
        
        comp[i,j] <- switch(stability.metric,
                            jaccard = {jaccard(features[,i], features[,j])},
                            sorenson = {sorenson(features[,i], features[,j])},
                            kuncheva = {kuncheva(features[,i], features[,j], nc)},
                            ochiai = {ochiai(features[,i], features[,j])},
                            pof = {pof(features[,i], features[,j])},
                            spearman = {spearman(features[,i], features[,j])},
                            canberra = {canberra(features[,i], features[,j])}
                            )
      }
    }
    
    # remove blank column (not sure why needed above, possibly remove)
    comp <- comp[,-1]
    
    # set row and column names
    colnames(comp) <- paste("Resample", 2:k, sep=".")
    rownames(comp) <- paste("Resample", 1:(k-1), sep=".")
    
    # Average all pairwise comparisons
    total <- round(2*sum(comp)/(k*(k-1)), 2)
    
    out <- list(comparisons = comp,
                overall = total)
    out  
  }

#' @title Pairwise Model Stability Metrics
#' @description Conducts all pairwise comparisons of each model's selected features selected 
#' following bootstrapping.  Also known as the function perturbation ensemble approach
#' @param features A matrix of selected features
#' @param stability.metric string indicating the type of stability metric.
#' Avialable options are \code{"jaccard"} (Jaccard Index/Tanimoto Distance),
#'  \code{"sorensen"} (Dice-Sorensen's Index), \code{"ochiai"} (Ochiai's Index),
#'  \code{"pof"} (Percent of Overlapping Features), \code{"kuncheva"} (Kuncheva's Stability Measures),
#'  \code{"spearman"} (Spearman Rank Correlation), and \code{"canberra"} (Canberra Distance)
#' @param m Number of models fit
#' @param k Number of resamples run
#' @return A list is returned containing:
#'  \item{comparisons}{Matrix of pairwise comparisons}
#'  \item{overall}{The average of all pairwise comparisons}
#' @author Charles Determan Jr
#' @references He. Z. & Weichuan Y. (2010) \emph{Stable feature selection for biomarker discovery}. 
#' Computational Biology and Chemistry 34 215-225.
pairwise.model.stability <- 
  function(features, stability.metric, m, k)
  {
    # column x of list one ~= to column x of list two
    # nro = number of pairwise comparisons
    nro = m*(m-1)/2
    tmp <- matrix(0, nrow=nro, ncol=k)
    tmp.names <- names(features)
    
    # used for identifying comparisons
    vec <- vector()
    for(i in 1:nro){
      vec[i] <- paste(tmp.names[i], ".vs.", tmp.names[i+1], sep = "")
    }
    
    a <- vector("list", k)
    b <- vector("list", k)
    
    # Retest this with fs.stability()
    for(i in 1:(m-1)){
      for(g in 1:k){
        #a[[g]] <- features[[i]]$features.selected[,g]
        #b[[g]] <- features[[i+1]]$features.selected[,g]
        a[[g]] <- features[[i]][,g]
        b[[g]] <- features[[i+1]][,g]
      }
    }
    
    tmp <- t(data.frame(mapply(a, FUN = stability.metric, y = b)))
    
    colnames(tmp) <- paste("Resample.", 1:k, sep = "")
    rownames(tmp) <- vec
    
    # Average all pairwise comparisons and resamples
    total <- round(sum(tmp)/k*nro, 2)
    out <- list(comparisons = tmp,
                overall = total)
    out  
  }

