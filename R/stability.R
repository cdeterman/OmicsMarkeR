
#####################################
### Stability/Similarity Measures ###
#####################################

#' @title Spearman Rank Correlation Coefficient
#' @description Calculates spearman rank correlation between two vectors
#' @param x numeric vector of ranks
#' @param y numeric vector of ranks with compatible length to x
#' @return Returns the spearman rank coefficient for the two vectors
#' @export

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
  
#' @title Canberra Distance
#' @description Calculates canberra distance between two vectors.  In brief, 
#' the higher the canberra distance the greater the 'distance' between the 
#' two vectors (i.e. they are less similar).
#' @param x numeric vector of ranks
#' @param y numeric vector of ranks with compatible length to x
#' @return Returns the canberra distance for the two vectors
#' @author Charles E. Determan Jr.
#' @references Jurman G., Merler S., Barla A., Paoli S., Galea A., & Furlanello C.
#' (2008) \emph{Algebraic stability indicators for ranked lists in molecular
#' profiling}. Bioinformatics 24(2): 258-264.
#' 
#' He. Z. & Weichuan Y. (2010) \emph{Stable feature selection for biomarker discovery}. 
#' Computational Biology and Chemistry 34 215-225.
#' @export

canberra <- function(x,y){
  x <- as.vector(as.numeric(as.character(x)))
  y <- as.vector(as.numeric(as.character(y)))
  out <- sum(abs(x-y)/(x+y))
  out
}

#' @title Jaccard Index
#' @description Calculates jaccard index between two vectors of features.  In brief, 
#' the closer to 1 the more similar the vectors.  The two vectors may have
#' an arbitrary cardinality (i.e. don't need same length).  Also known as 
#' the Tanimoto distance metric.  Defined as the size of the vectors' intersection
#' divided by the size of the union of the vectors.
#' @param x vector of feature names
#' @param y vector of feature names
#' @return Returns the jaccard index for the two vectors. It takes values in [0,1], 
#' with 0 meaning no overlap between two sets and 1 meaning two sets are identical.
#' @author Charles E. Determan Jr.
#' @references Jaccard P. (1908) \emph{Nouvelles recherches sur la distribution florale}. 
#' Bull. Soc. Vaudoise Sci. Nat. 44: 223-270.
#' 
#' Real R. & Vargas J.M. (1996) \emph{The Probabilistic Basis of Jaccard's Index of Similarity}
#' Systematic Biology 45(3): 380-385.
#' 
#' He. Z. & Weichuan Y. (2010) \emph{Stable feature selection for biomarker discovery}. 
#' Computational Biology and Chemistry 34 215-225.
#' @seealso \code{\link{kuncheva}}, \code{\link{sorensen}}, \code{\link{ochiai}},
#' \code{\link{pof}}, \code{\link{pairwise.stability}}, \code{\link{pairwise.model.stability}}
#' @export
jaccard <- function(x,y){
  x <- as.vector(x)
  y <- as.vector(y)
  index <- length(intersect(x,y))/length(union(x, y))
  return(index)
}

#' @title Dice-Sorensen's Index
#' @description Calculates Dice-Sorensen's index between two vectors of features.  In brief, 
#' the closer to 1 the more similar the vectors.  The two vectors may have
#' an arbitrary cardinality (i.e. don't need same length).  Very similar to the
#' Jaccard Index \code{\link{jaccard}} but Dice-Sorensen is the harmonic mean of
#' the ratio.
#' @param x vector of feature names
#' @param y vector of feature names
#' @return Returns the Dice-Sorensen's Index for the two vectors. It takes values in [0,1], with 0 meaning no overlap 
#' between two sets and 1 meaning two sets are identical.
#' @author Charles E. Determan Jr.
#' @references Sorensen T. (1957) \emph{A method of establishing roups of equal amplitude
#' in plant sociology based on similarity of species and its application to analyses of
#' the vegetation on Danish commons}. 
#' Kongelige Danske Videnskabernes Selskab. 5(4): 1-34.
#' 
#' Dice, Lee R. (1945) \emph{Measures of the Amount of Ecologic Association
#' Between Species}. Ecology 26 (3): 297-302. doi:10.2307/1932409
#' 
#' He. Z. & Weichuan Y. (2010) \emph{Stable feature selection for biomarker discovery}. 
#' Computational Biology and Chemistry 34 215-225.
#' @seealso \code{\link{kuncheva}}, \code{\link{sorensen}}, \code{\link{ochiai}},
#' \code{\link{pof}}, \code{\link{pairwise.stability}}, \code{\link{pairwise.model.stability}}
#' @export

sorensen <- function(x,y){
  x <- as.vector(x)
  y <- as.vector(y)
  index <- 2*(length(intersect(x,y)))/(2*(length(intersect(x,y)))+length(setdiff(x,y))+length(setdiff(y,x)))
  return(index)
}

#' @title Ochiai's Index
#' @description Calculates Ochiai's index between two vectors of features.  In brief, 
#' the closer to 1 the more similar the vectors.  The two vectors may have
#' an arbitrary cardinality (i.e. don't need same length).  Very similar to the
#' Jaccard Index \code{\link{jaccard}} but Ochiai is a geometric means of the
#' ratio.  
#' @param x vector of feature names
#' @param y vector of feature names
#' @return Returns the Ochiai Index for the two vectors. It takes values in [0,1], with 0 meaning no overlap 
#' between two sets and 1 meaning two sets are identical.
#' @author Charles E. Determan Jr.
#' @references Ochiai A. (1957) \emph{Zoogeographical studies on the soleoid 
#' fishes found in Japan and its neigbouring regions}. 
#' Bulletin of the Japanese Society of Scientific Fisheries. 22: 526-530.
#' 
#' Zucknick M., Richardson S., & Stronach E.A. (2008) \emph{Comparing the characteristics of gene expression 
#' profiles derived by univariate and multivariate classification methods}. 
#' Statistical Applications in Genetics and Molecular Biology. 7(1): Article 7. doi:10.2202/1544-6115.1307
#' 
#' He. Z. & Weichuan Y. (2010) \emph{Stable feature selection for biomarker discovery}. 
#' Computational Biology and Chemistry 34 215-225.
#' @seealso \code{\link{kuncheva}}, \code{\link{sorensen}}, \code{\link{ochiai}},
#' \code{\link{pof}}, \code{\link{pairwise.stability}}, \code{\link{pairwise.model.stability}}
#' @export

ochiai <- function(x,y){
  x <- as.vector(x)
  y <- as.vector(y)
  index <- 2*(length(intersect(x,y)))/(sqrt(length(x)*length(y)))
  return(index)
}

#' @title Percentage of Overlapping Features
#' @description Calculates percent of overlapping features between two vectors of features.  In brief, 
#' the closer to 1 the more similar the vectors.  The two vectors may have
#' an arbitrary cardinality (i.e. don't need same length).
#' @param x vector of feature names
#' @param y vector of feature names
#' @return Returns the percent of overlapping features for the two vectors. It takes values in [0,1], with 0 meaning no overlap 
#' between two sets and 1 meaning two sets are identical.
#' @author Charles E. Determan Jr.
#' @references Shi L., et al. (2005) \emph{Cross-platform comparability of microarray
#' technology: intra-platform consistency and appropriate data analysis procedures are essential}. 
#' BMC Bioinformatics. 6 (Suppl. 2) S12.

#' He. Z. & Weichuan Y. (2010) \emph{Stable feature selection for biomarker discovery}. 
#' Computational Biology and Chemistry 34 215-225.
#' @seealso \code{\link{kuncheva}}, \code{\link{sorensen}}, \code{\link{ochiai}},
#' \code{\link{pof}}, \code{\link{pairwise.stability}}, \code{\link{pairwise.model.stability}}
#' @export

pof <- function(x,y){
  x <- as.vector(x)
  y <- as.vector(y)
  index <- length(intersect(x,y))/length(x)
  return(index)
}

#' @title Kuncheva's Index
#' @description Calculates Kuncheva's index between two vectors of features.  In brief, 
#' the closer to 1 the more similar the vectors.  The two vectors must have the same 
#' cardinality (i.e. same length).
#' @param x vector of feature names
#' @param y vector of feature names
#' @param num.features total number of features in the original dataset
#' @return Returns the Kuncheva Index for the two vectors. It takes values in [0,1], with 0 meaning no overlap 
#' between two sets and 1 meaning two sets are identical.
#' @author Charles E. Determan Jr.
#' @references Kuncheva L. (2007) \emph{A stability index for feature selection}. 
#' Proceedings of the 25th IASTED International Multi-Conference: Artificial Intelligence and 
#' Applications. pp. 390-395.
#' 
#' He. Z. & Weichuan Y. (2010) \emph{Stable feature selection for biomarker discovery}. 
#' Computational Biology and Chemistry 34 215-225.
#' @seealso \code{\link{kuncheva}}, \code{\link{sorensen}}, \code{\link{ochiai}},
#' \code{\link{pof}}, \code{\link{pairwise.stability}}, \code{\link{pairwise.model.stability}}
#' @export

kuncheva <- function(x,             # list of features in first run
                     y,             # list of features in second run
                     num.features   # total number of features in entire dataset
                     )
  {
  # m/N = total number of features
  # s/c = length of features - assumes equal length of x and y
  # r = number of common elements in both signatures
  # (r*N - s^2)/s*(N-s)
  
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
#'  @param nc Number of variables in original dataset
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
                            sorensen = {sorensen(features[,i], features[,j])},
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
    for(i in 1:(m-1)){
      for(g in 1:k){
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

