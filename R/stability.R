
#####################################
### Stability/Similarity Measures ###
#####################################

#' @title Spearman Rank Correlation Coefficient
#' @description Calculates spearman rank correlation between two vectors
#' @param x numeric vector of ranks
#' @param y numeric vector of ranks with compatible length to x
#' @return Returns the spearman rank coefficient for the two vectors
#' @example inst/examples/spearman.R
#' @export

spearman <- function(x,y){
    assert_is_numeric(x)
    assert_is_numeric(y)
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
#' @note The \code{\link{canberra_stability}} function is used
#' internally to return the canberra metric.
#' @author Charles E. Determan Jr.
#' @references Jurman G., Merler S., Barla A., Paoli S., Galea A., & 
#' Furlanello C. (2008) \emph{Algebraic stability indicators for ranked lists 
#' in molecular profiling}. Bioinformatics 24(2): 258-264.
#' 
#' He. Z. & Weichuan Y. (2010) \emph{Stable feature selection for biomarker 
#' discovery}. Computational Biology and Chemistry 34 215-225.
#' @example inst/examples/canberra.R
#' @export

canberra <- function(x,y){
    assert_is_numeric(x)
    assert_is_numeric(y)
    if(length(x) != length(y)){
        stop("\n Error: feature lengths must be same for Canberra Distance")
    }
    out <- sum(abs(x-y)/(x+y))
    out
}

expected_canberra <- function(p, k=p){
    out <- (((k+1)*(2*p - k)*log(4))/p) - ((2*k*p + 3*p - k - k^2)/p)
    
    #out <- (log(4) - 1)*p + log(4) - 2
    return(out)
}


#' @title Canberra Stability
#' @description Calculates canberra stability between two ranked lists.  In 
#' brief, the raw canberra distance is scaled to a [0,1] distribution by the
#' maximum canberra metric.  Lastly, this value is subtracted from 1 to provide 
#' the same interpretation as the other stability metrics whereby 1 is 
#' identical and 0 is no stability.
#' @param x numeric vector of ranks
#' @param y numeric vector of ranks with compatible length to x
#' @return Returns the canberra stability for the two vectors
#' @author Charles E. Determan Jr.
#' @references Jurman G., Merler S., Barla A., Paoli S., Galea A., & 
#' Furlanello C. (2008) \emph{Algebraic stability indicators for ranked lists 
#' in molecular profiling}. Bioinformatics 24(2): 258-264.
#' 
#' He. Z. & Weichuan Y. (2010) \emph{Stable feature selection for biomarker 
#' discovery}. Computational Biology and Chemistry 34 215-225.
#' @example inst/examples/canberra.R
#' @export
canberra_stability <- function(x,y){
    d <- sum(abs(x-y)/(x+y))
    # scale to [0,1] range
    # (x - xmin)/(xmax - xmin)
    d_max <- canberra(seq(length(x)), rev(seq(length(x))))
    out <- 1 - (d - 0)/(d_max - 0)
    return(out)
}

#' @title Jaccard Index
#' @description Calculates jaccard index between two vectors of features.  
#' In brief, the closer to 1 the more similar the vectors.  The two vectors 
#' may have an arbitrary cardinality (i.e. don't need same length).  Also 
#' known as the Tanimoto distance metric.  Defined as the size of the vectors' 
#' intersection divided by the size of the union of the vectors.
#' @param x vector of feature names
#' @param y vector of feature names
#' @return Returns the jaccard index for the two vectors. It takes values 
#' in [0,1], with 0 meaning no overlap between two sets and 1 meaning two sets 
#' are identical.
#' @author Charles E. Determan Jr.
#' @references Jaccard P. (1908) \emph{Nouvelles recherches sur la 
#' distribution florale}. Bull. Soc. Vaudoise Sci. Nat. 44: 223-270.
#' 
#' Real R. & Vargas J.M. (1996) \emph{The Probabilistic Basis of Jaccard's 
#' Index of Similarity} Systematic Biology 45(3): 380-385.
#' 
#' He. Z. & Weichuan Y. (2010) \emph{Stable feature selection for 
#' biomarker discovery}. Computational Biology and Chemistry 34 215-225.
#' @seealso \code{\link{kuncheva}}, \code{\link{sorensen}}, 
#' \code{\link{ochiai}}, \code{\link{pof}}, \code{\link{pairwise.stability}}, 
#' \code{\link{pairwise.model.stability}}
#' @example inst/examples/jaccard.R
#' @export
jaccard <- function(x,y){
    assert_is_character(x)
    assert_is_character(y)
    
    x <- as.vector(x)
    y <- as.vector(y)
    index <- length(intersect(x,y))/length(union(x, y))
    return(index)
}

#' @title Dice-Sorensen's Index
#' @description Calculates Dice-Sorensen's index between two vectors of 
#' features.  In brief, the closer to 1 the more similar the vectors.  
#' The two vectors may have an arbitrary cardinality (i.e. don't need 
#' same length).  Very similar to the Jaccard Index \code{\link{jaccard}} 
#' but Dice-Sorensen is the harmonic mean of the ratio.
#' @param x vector of feature names
#' @param y vector of feature names
#' @return Returns the Dice-Sorensen's Index for the two vectors. It takes 
#' values in [0,1], with 0 meaning no overlap between two sets and 1 meaning 
#' two sets are identical.
#' @author Charles E. Determan Jr.
#' @references Sorensen T. (1948) \emph{A method of establishing roups of 
#' equal amplitude in plant sociology based on similarity of species and 
#' its application to analyses of the vegetation on Danish commons}. 
#' Kongelige Danske Videnskabernes Selskab. 5(4): 1-34.
#' 
#' Dice, Lee R. (1945) \emph{Measures of the Amount of Ecologic Association
#' Between Species}. Ecology 26 (3): 297-302. doi:10.2307/1932409
#' 
#' He. Z. & Weichuan Y. (2010) \emph{Stable feature selection for 
#' biomarker discovery}. Computational Biology and Chemistry 34 215-225.
#' @seealso \code{\link{kuncheva}}, \code{\link{sorensen}}, 
#' \code{\link{ochiai}}, \code{\link{pof}}, \code{\link{pairwise.stability}}, 
#' \code{\link{pairwise.model.stability}}
#' @example inst/examples/sorensen.R
#' @export

sorensen <- function(x,y){
    assert_is_character(x)
    assert_is_character(y)
    index <- 
        2*(length(intersect(x,y)))/(2*(length(intersect(x,y)))+
                                        length(setdiff(x,y))+
                                        length(setdiff(y,x)))
    return(index)
}

#' @title Ochiai's Index
#' @description Calculates Ochiai's index between two vectors of features.  
#' In brief, the closer to 1 the more similar the vectors.  The two vectors 
#' may have an arbitrary cardinality (i.e. don't need same length).  Very 
#' similar to the Jaccard Index \code{\link{jaccard}} but Ochiai is a 
#' geometric means of the ratio.  
#' @param x Character vector of feature names
#' @param y Character vector of feature names
#' @return Returns the Ochiai Index for the two vectors. It takes values 
#' in [0,1], with 0 meaning no overlap between two sets and 1 meaning two 
#' sets are identical.
#' @author Charles E. Determan Jr.
#' @references Ochiai A. (1957) \emph{Zoogeographical studies on the soleoid 
#' fishes found in Japan and its neigbouring regions}. 
#' Bulletin of the Japanese Society of Scientific Fisheries. 22: 526-530.
#' 
#' Zucknick M., Richardson S., & Stronach E.A. (2008) \emph{Comparing the 
#' characteristics of gene expression profiles derived by univariate and 
#' multivariate classification methods}. Statistical Applications in Genetics 
#' and Molecular Biology. 7(1): Article 7. doi:10.2202/1544-6115.1307
#' 
#' He. Z. & Weichuan Y. (2010) \emph{Stable feature selection for 
#' biomarker discovery}. Computational Biology and Chemistry 34 215-225.
#' @seealso \code{\link{kuncheva}}, \code{\link{sorensen}}, 
#' \code{\link{ochiai}}, \code{\link{pof}}, \code{\link{pairwise.stability}}, 
#' \code{\link{pairwise.model.stability}}
#' @example inst/examples/ochiai.R
#' @export

ochiai <- function(x,y){
    assert_is_character(x)
    assert_is_character(y)
    index <- length(intersect(x,y))/(sqrt(length(x)*length(y)))
    return(index)
}

#' @title Percentage of Overlapping Features
#' @description Calculates percent of overlapping features between two vectors 
#' of features.  In brief, the closer to 1 the more similar the vectors.  
#' The two vectors may have an arbitrary cardinality (i.e. don't need 
#' same length).
#' @param x Character vector of feature names
#' @param y Character vector of feature names
#' @return Returns the percent of overlapping features for the two vectors. 
#' It takes values in [0,1], with 0 meaning no overlap between two sets and 1 
#' meaning two sets are identical.
#' @author Charles E. Determan Jr.
#' @references Shi L., et al. (2005) \emph{Cross-platform comparability 
#' of microarray technology: intra-platform consistency and appropriate data 
#' analysis procedures are essential}. BMC Bioinformatics. 6 (Suppl. 2) S12.

#' He. Z. & Weichuan Y. (2010) \emph{Stable feature selection for 
#' biomarker discovery}. Computational Biology and Chemistry 34 215-225.
#' @seealso \code{\link{kuncheva}}, \code{\link{sorensen}}, 
#' \code{\link{ochiai}}, \code{\link{pof}}, \code{\link{pairwise.stability}}, 
#' \code{\link{pairwise.model.stability}}
#' @example inst/examples/pof.R
#' @export

pof <- function(x,y){
    assert_is_character(x)
    assert_is_character(y)
    index <- length(intersect(x,y))/length(x)
    return(index)
}

#' @title Kuncheva's Index
#' @description Calculates Kuncheva's index between two vectors of features.  
#' In brief, the closer to 1 the more similar the vectors.  The two vectors 
#' must have the same cardinality (i.e. same length).
#' @param x Character vector of feature names
#' @param y Character vector of feature names
#' @param num.features total number of features in the original dataset
#' @return Returns the Kuncheva Index for the two vectors. It takes values 
#' in [0,1], with 0 meaning no overlap between two sets and 1 meaning two 
#' sets are identical.
#' @note The returned Kuncheva Index has been scaled from its original [-1,1]
#' range to [0,1] in order to make it compatible with RPT.
#' @author Charles E. Determan Jr.
#' @references Kuncheva L. (2007) \emph{A stability index for feature 
#' selection}. Proceedings of the 25th IASTED International Multi-Conference: 
#' Artificial Intelligence and Applications. pp. 390-395.
#' 
#' He. Z. & Weichuan Y. (2010) \emph{Stable feature selection for 
#' biomarker discovery}. Computational Biology and Chemistry 34 215-225.
#' @seealso \code{\link{kuncheva}}, \code{\link{sorensen}}, 
#' \code{\link{ochiai}}, \code{\link{pof}}, \code{\link{pairwise.stability}}, 
#' \code{\link{pairwise.model.stability}}
#' @example inst/examples/kuncheva.R
#' @export

kuncheva <- function(x,             
                     y, 
                     num.features 
)
{
    # m/N = total number of features
    # s/c = length of features - assumes equal length of x and y
    # r = number of common elements in both signatures
    # (r*N - s^2)/s*(N-s)
    
    assert_is_not_null(num.features)
    
    assert_is_character(x)
    assert_is_character(y)
    if(length(x) != length(y)){
        stop("\n Error: Kuncheva requires same number of features in subset.
             \nYou must specifiy number of features to be selected")    
    }
    
    r <- length(intersect(x,y))
    d <- num.features
    k <- length(x)      # cardinality
    
    index <- (r - ((k^2)/d))/(k - ((k^2)/d))
    
    # scale to [0,1] range
    # (x - xmin)/(xmax - xmin)
    index <- (index - -1)/(1 - -1)
    # scale index to 0,1 range
    return(index)
    }

#' @title Pairwise Stability Metrics
#' @description Conducts all pairwise comparisons of features selected 
#' following bootstrapping. Also known as the data perturbation ensemble 
#' approach.
#' @param features A matrix of selected features
#' @param stability.metric string indicating the type of stability metric.
#' Available options are \code{"jaccard"} (Jaccard Index/Tanimoto Distance),
#'  \code{"sorensen"} (Dice-Sorensen's Index), \code{"ochiai"} (Ochiai's Index),
#'  \code{"pof"} (Percent of Overlapping Features), \code{"kuncheva"} 
#'  (Kuncheva's Stability Measures), \code{"spearman"} (Spearman Rank 
#'  Correlation), and \code{"canberra"} (Canberra Distance)
#'  @param nc Number of variables in original dataset
#' @return A list is returned containing:
#' \item{comparisons}{Matrix of pairwise comparisons}
#' \item{overall}{The average of all pairwise comparisons}
#' @author Charles Determan Jr
#' @references He. Z. & Weichuan Y. (2010) \emph{Stable feature selection for 
#' biomarker discovery}. Computational Biology and Chemistry 34 215-225.
#' @example inst/examples/pairwise.stability.R
#' @export

pairwise.stability <- 
    function(features, 
             stability.metric, 
             nc
    )
    {
        assert_is_matrix(features)
        assert_is_character(stability.metric)
        
        if(stability.metric == "kuncheva"){
            assert_is_not_null(c(nc))
            assert_is_numeric(nc)
        }
        
        k <- ncol(features)
        comp <- matrix(0, nrow=k-1, ncol=k)
        for (i in 1:(k-1)){
            for (j in (i+1):k){
                
                comp[i,j] <- switch(
                    stability.metric,
                    jaccard = {jaccard(features[,i], features[,j])},
                    sorensen = {sorensen(features[,i], features[,j])},
                    kuncheva = {kuncheva(features[,i], features[,j], nc)},
                    ochiai = {ochiai(features[,i], features[,j])},
                    pof = {pof(features[,i], features[,j])},
                    spearman = {spearman(features[,i], features[,j])},
                    canberra = {canberra_stability(features[,i], features[,j])}
                )
            }
        }
        
        # remove blank column (not sure why needed above, possibly remove)
        comp <- comp[,-1]
        
        # set row and column names
        if(k > 2){
            colnames(comp) <- paste("Resample", 2:k, sep=".")
            rownames(comp) <- paste("Resample", 1:(k-1), sep=".")
        }
        
        # Average all pairwise comparisons
        total <- round(2*sum(comp)/(k*(k-1)), 2)
        
        out <- list(comparisons = comp,
                    overall = total)
        out  
    }

#' @title Pairwise Model Stability Metrics
#' @description Conducts all pairwise comparisons of each model's selected 
#' features selected  following bootstrapping.  Also known as the function 
#' perturbation ensemble approach
#' @param features A matrix of selected features
#' @param stability.metric string indicating the type of stability metric.
#' Avialable options are \code{"jaccard"} (Jaccard Index/Tanimoto Distance),
#'  \code{"sorensen"} (Dice-Sorensen's Index), \code{"ochiai"} (Ochiai's Index),
#'  \code{"pof"} (Percent of Overlapping Features), \code{"kuncheva"} 
#'  (Kuncheva's Stability Measures), \code{"spearman"} (Spearman Rank 
#'  Correlation), and \code{"canberra"} (Canberra Distance)
#' @param nc Number of original features
#' @return A list is returned containing:
#'  \item{comparisons}{Matrix of pairwise comparisons}
#'  \item{overall}{The average of all pairwise comparisons}
#' @author Charles Determan Jr
#' @references He. Z. & Weichuan Y. (2010) \emph{Stable feature selection for 
#' biomarker discovery}. Computational Biology and Chemistry 34 215-225.
#' @seealso \code{\link{pairwise.stability}}
#' @example inst/examples/pairwise.model.stability.R
#' @import plyr
#' @export

pairwise.model.stability <- 
    function(features, stability.metric, nc)
    {
        if(any(!sapply(features, is.data.frame))){
            features <- lapply(features, as.data.frame)
        }
        # column x of list one ~= to column x of list two
        # nro = number of pairwise comparisons
        m <- length(features)
        k <- unique(unlist(lapply(features, ncol)))
        if(length(k) > 1){
            stop("\n Error: The number of resamples must be equal")
        }
        nro <- m*(m-1)/2
        tmp <- matrix(0, nrow=nro, ncol=k)
        tmp.names <- names(features)
        
        # used for identifying comparisons
        vec <- matrix(0, nrow=m-1, ncol=m)
        #vec <- vector()
        for(i in 1:(m-1)){
            for(j in (i+1):m){
                vec[i,j] <- paste(tmp.names[i], ".vs.", tmp.names[j], sep = "")
            }
        }
        vec.comps <- as.vector(t(vec))[as.vector(t(vec)) != 0]    
        
        if(length(unique(unlist(lapply(features, 
                                       function(x) dim(x)[1])))) == 1){
            model.features <- vector("list", k)
            for(i in seq(k)){
                model.features[[i]] <- sapply(features, `[[`, i)
            }
        }else{
            model.features <- vector("list", k)
            for(c in seq(k)){
                tmp <- 
                    lapply(features, 
                           FUN = function(x){
                               as.data.frame(t(data.frame(x[,c])))
                               })
                model.features[[c]] <- t(rbind.fill(tmp))
            } 
        }
        
        extract.pairs <- function(tmp, m){
            s <- seq(m-1)
            out <- c()
            for(i in s){
                tmp.out <- as.matrix(tmp)[i,i:max(s)]
                out <- c(out, tmp.out)
            }
            out
        }
        
        tmp <- 
            lapply(model.features, 
                   FUN = function(x){
                       pairwise.stability(x, stability.metric, nc)$comparisons
                   })
        tmp.dat <- sapply(tmp, FUN = function(x,m) extract.pairs(x, m), m = m)
        
        if(!is.matrix(tmp.dat)){
            tmp.dat <- t(as.matrix(tmp.dat))
        }
        
        colnames(tmp.dat) <- paste("Resample.", 1:k, sep = "")
        rownames(tmp.dat) <- vec.comps
        
        # Average all pairwise comparisons and resamples
        total <- round(sum(tmp.dat)/(k*nro), 2)
        out <- list(comparisons = tmp.dat,
                    overall = total)
        out  
    }

#' @title Robustness-Performance Trade-Off
#' @description A variation on the F-measure (precision and recall) to
#' assess robustness versus classification performance.
#' @param stability Stability metric i.e. result from jaccard, sorensen, etc.
#' @param performance Model performance e.g. accuracy
#' @param beta Relative of importance of stability versus performance.
#' Default \code{beta = 1} treats stability and performance equally.
#' @return Harmonic mean of robustness and classification performance
#' @references Saeys Y., Abeel T., et. al. (2008) \emph{Machine Learning and 
#' Knowledge Discovery in Databases}. 313-325. 
#' http://link.springer.com/chapter/10.1007/978-3-540-87481-2_21
#' @example inst/examples/RPT.R
#' @export

RPT <- 
    function(stability, performance, beta = 1){
        assert_is_numeric(stability)
        assert_is_numeric(performance)
        assert_is_numeric(beta)
        assert_is_in_closed_range(stability, lower=0, upper=1)
        assert_is_in_closed_range(performance, lower=0, upper=1)
        assert_is_positive(beta)
        
        num <- ((beta^2)+1)*stability*performance
        denom <- (beta^2)*stability + performance
        rpt <- num/denom
        rpt
    }
