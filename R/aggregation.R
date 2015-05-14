###########################
### Aggregation Methods ###
###########################

## Notation
# s = selected number of top ranking features
# r = feature rank
# f = matrix of ranked features by bootstrap

#' @title Complete Linear Aggregation
#' @description Compiles matrix of ranked features via complete 
#' linear aggregation
#' @param efs A matrix of selected features
#' @param f The number of features desired.  If rank correlation 
#' desired, f = NULL
#' @return \item{agg}{Aggregated list of features}
#' @seealso \code{\link{ES}}, \code{\link{EM}}, \code{\link{EE}}
#' @author Charles Determan Jr
#' @references Abeel T., Helleputte T., Van de Peer Y., Dupont P., Saeys Y. 
#' (2010) \emph{Robust biomarker identification for cancer diagnosis with 
#' ensemble feature selection methods}. Bioinformatics 26:3 392-398.
# ' @export

CLA <- function(efs, f){
    
    assert_is_matrix(efs)
    
    # sum ranks over all bootstraps (i.e. colsums)
    agg <- rank(rowSums(efs))
    if(!is.null(f)){
        agg <- head(sort(agg), f)
    }else{
        agg <- sort(agg)
    }
    ## return features with ranks
    as.matrix(agg)
}

#' @title Ensemble Mean Aggregation
#' @description Compiles matrix of ranked features via ensemble 
#' mean aggregation
#' @param efs A matrix of selected features
#' @param f The number of features desired.  If rank correlation desired,
#' f = NULL
#' @return \item{agg}{Aggregated list of features}
#' @seealso \code{\link{CLA}}, \code{\link{ES}}, \code{\link{EE}}
#' @author Charles Determan Jr
#' @references Abeel T., Helleputte T., Van de Peer Y., Dupont P., Saeys Y. 
#' (2010) \emph{Robust biomarker identification for cancer diagnosis with 
#' ensemble feature selection methods}. Bioinformatics 26:3 392-398.
# ' @export

EM <- function(efs, f){
    assert_is_matrix(efs)
    
    # average rank over the bootstraps
    agg <- rank(rowMeans(efs, na.rm = FALSE))
    
    # return list of 's' features with highest scores
    if(!is.null(f)){
        agg <- head(sort(agg), f)
    }else{
        agg <- sort(agg)
    }
    
    ## return features with ranks
    as.matrix(agg)
}

#' @title Ensemble Stability Aggregation
#' @description Compiles matrix of ranked features via ensemble 
#' stability aggregation
#' @param efs A matrix of selected features
#' @param f The number of features desired.  If rank correlation desired, 
#' f = NULL
#' @return \item{agg}{Aggregated list of features}
#' @seealso \code{\link{CLA}}, \code{\link{EM}}, \code{\link{EE}}
#' @author Charles Determan Jr
#' @references Meinshausen N., Buhlmann P. (2010) \emph{Stability selection}. 
#' J.R. Statist. Soc. B. 72:4 417-473.
# ' @export

ES <- function(efs, f){
    assert_is_matrix(efs)
    
    # measures percentage of bootstrap samples for which the gene ranks 
    # in the top 's'
    # i.e. f(r) = 1 if r <= s, otherwise = 0
    
    agg <- rank(-apply(efs, 1, FUN = function(x){
        length(which(x <= 10)==TRUE)/length(x)})
    )
    # return list of 's' features with highest scores
    
    if(!is.null(f)){
        agg <- head(sort(agg), f)
    }else{
        agg <- sort(agg)
    }
    
    as.matrix(agg)
}

#' @title Ensemble Exponential Aggregation
#' @description Compiles matrix of ranked features via ensemble 
#' exponential aggregation
#' @param efs A matrix of selected features
#' @param f The number of features desired.  If rank correlation desired, 
#' f = NULL
#' @return \item{agg}{Aggregated list of features}
#' @seealso \code{\link{CLA}}, \code{\link{ES}}, \code{\link{EM}},
#' @author Charles Determan Jr
#' @references Haury A., Gestraud P., Vert J. (2011) \emph{The Influence of 
#' Features Selection Methods on Accuracy, Stability, and Interpretability 
#' of Molecular Signatures}.  PLoS ONE 6(12) e28210. doi: 
#' 10.1371/journal.pone.0028210
# ' @export

EE <- function(efs, f){
    assert_is_matrix(efs)
    assert_is_numeric(f)
    
    # average an exponentially decreasing function of the rank
    # f(r) = exp(-r/s)
    # return list of 's' features with highest scores
    
    agg <- rank(-rowMeans(exp(-efs/f)))
    agg <- head(sort(agg), f)
    as.matrix(agg)
    }

#' @title Feature Aggregation
#' @description Compiles matrix of ranked features via user defined 'metric' 
#' @param efs A matrix of selected features
#' @param metric string indicating the type of aggregation.
#' Avialable options are \code{"CLA"} (Complete Linear),
#'  \code{"EM"} (Ensemble Mean), \code{"ES"} (Ensemble Stability), and
#'  \code{"EE"} (Ensemble Exponential)
#' @param f The number of features desired.  Default \code{f = NULL}
#' @return \item{agg}{Aggregated list of features}
#' @seealso \code{\link{CLA}}, \code{\link{ES}}, \code{\link{EM}},
#' \code{\link{EE}}
#' @author Charles Determan Jr
#' @references Abeel T., Helleputte T., Van de Peer Y., Dupont P., Saeys Y. 
#' (2010) \emph{Robust biomarker identification for cancer diagnosis with 
#' ensemble feature selection methods}. Bioinformatics 26(3) 392-398.
#' 
#' Meinshausen N., Buhlmann P. (2010) \emph{Stability selection}. 
#' J.R. Statist. Soc. B. 72(4) 417-473.
#' 
#' Haury A., Gestraud P., Vert J. (2011) \emph{The Influence of Features 
#' Selection Methods on Accuracy, Stability, and Interpretability of 
#' Molecular Signatures}.  PLoS ONE 6(12) e28210. doi: 
#' 10.1371/journal.pone.0028210.
#' @example /inst/examples/aggregation.R
#' @export

aggregation <- function(efs,
                        metric,
                        f = NULL
)
{
    assert_is_matrix(efs)
    
    agg <- switch(metric,
                  CLA = {CLA(efs, f)},
                  EM = {EM(efs, f)},
                  ES = {ES(efs, f)},
                  EE = {EE(efs, f)}
    )
    agg
}
