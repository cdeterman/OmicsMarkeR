#' @title Feature Extraction
#' @description Extracts features from models that have been previously fit.
#' @param x Previously fitted model
#' @param dat Numeric variable data used for fitted models 
#' (In appropriate format)
#' @param grp Vector of training classes
#' @param method String indicating the INDIVIDUAL model being extracted from
#' @param model.features Logical argument dictating if features selected 
#' determined by models instead of user determined number of features.
#' @param bestTune If \code{model.features = TRUE}, must provide the parameter 
#' at which to extract features from the model.
#' @param f Number of features to subset
#' @param comp.catch An internal check for plsda models.  If the optimal model 
#' contains only 1 component, the ncomp paramter must be set to 2 for the 
#' model.  However, features are still extracted only from the first component.
#' @return Returns list of the features selected from the fitted model.
#' @import randomForest
#' @import e1071
#' @import gbm
#' @import pamr
#' @import glmnet
# ' @export

extract.features <-
    function(x, 
             dat = NULL,  
             grp = NULL,  
             method,    
             model.features = FALSE, 
             bestTune = NULL, 
             f, 
             comp.catch = NULL 
    )
    {
        features <- 
            switch(method,
                   plsda =
{
    if(!is.null(comp.catch)){
        l <- 1
    }else{
        l <- ncol(x[[1]]$VIP)  
    }
    
    if(model.features){
        # generally a VIP >= 1 is significant in the literature
        sub <- lapply(x, 
                      FUN = function(x) x$VIP[which(x$VIP[,l] >= 1),])
        mod.features <- 
            lapply(sub, 
                   FUN = function(x) names(sort(x[,l], decreasing= TRUE)))
    }else{
        if(is.null(f)){
            # plsda VIP = higher is more important so must 
            # reverse rank function with '-'
            mod.features <- lapply(x, FUN = function(x) rank(-x$VIP[,l]))
        }else{
            mod.features <- 
                lapply(x, 
                       FUN = function(x){
                           names(sort(x$VIP[,l], decreasing = TRUE))
                       })
            mod.features <- lapply(mod.features, head, n = f)
        }
    }
    
    ## Make a matrix from the list, with shorter 
    # vectors filled out with ""
    n <- max(sapply(mod.features, FUN = function(x) length(x)))
    ll <- lapply(mod.features, function(x) {
        c(as.character(x), rep("", times = n - length(x)))
    })
    out <- as.data.frame(do.call(cbind, ll)) 
    
    # collect ranked features
    list(features.selected = out)
},

gbm =  
{
    if(model.features){
        # generally a VIP >= 1 is significant in the literature
        vips <- 
            lapply(x, 
                   FUN = function(x){
                       relative.influence(x, n.trees = bestTune$.n.trees)
                   })
        sub <- lapply(vips, FUN = function(x) x[which(x >= 1)])
        mod.features <- 
            lapply(sub, FUN = function(x) names(sort(x, decreasing= TRUE)))
    }else{
        if(is.null(f)){
            # plsda VIP = higher is more important so must 
            # reverse rank function with '-'
            mod.features <- 
                lapply(x, 
                       FUN = function(x){
                           rank(-relative.influence(x, n.trees = x$n.trees))
                       }
                )
        }else{
            mod.features <- 
                lapply(x, 
                       FUN = function(x){
                           names(
                               relative.influence(x, 
                                                  n.trees = x$n.trees, 
                                                  sort.= TRUE))
                       })
            mod.features <- lapply(mod.features, head, n = f)
        }
    }
    
    ## Make a matrix from the list, with shorter vectors 
    # filled out with ""
    n <- max(sapply(mod.features, FUN = function(x) length(x)))
    ll <- lapply(mod.features, function(x) {
        c(as.character(x), rep("", times = n - length(x)))
    })
    out <- as.data.frame(do.call(cbind, ll)) 
    
    # collect ranked features
    list(features.selected = out)
},

rf =
{
    # Mean Decrease in Accuracy metric (type=1)
    # Gini Index (type=2)
    
    if(model.features){
        # generally a VIP >= 1 is significant in the literature
        vips <- lapply(x, FUN = function(x) importance(x, type = 1))
        sub <- 
            lapply(vips, FUN = function(x) x[which(x >= 1),,drop = FALSE])
        mod.features <- 
            lapply(sub, 
                   FUN = function(x){
                       names(
                           x[sort(x, 
                                  decreasing = TRUE, 
                                  index.return= TRUE)$ix,])
                   })
    }else{
        if(is.null(f)){
            # random forest VIP = higher is more important so 
            # must reverse rank function with '-'
            mod.features <- 
                lapply(x, FUN = function(x) rank(-importance(x, type = 1)))
        }else{
            vips <- lapply(x, FUN = function(x) importance(x, type = 1))
            mod.features <- 
                lapply(vips, 
                       FUN = function(x){
                           names(
                               x[sort(x, 
                                      decreasing = TRUE, 
                                      index.return= TRUE)$ix,]
                           )
                       })
            mod.features <- lapply(mod.features, head, n = f)
        }
    }
    
    ## Make a matrix from the list, with shorter vectors 
    # filled out with ""
    n <- max(sapply(mod.features, FUN = function(x) length(x)))
    ll <- lapply(mod.features, function(x) {
        c(as.character(x), rep("", times = n - length(x)))
    })
    out <- as.data.frame(do.call(cbind, ll)) 
    
    # collect ranked features
    list(features.selected = out)
},                  

svm =   
{
    best.C <- c(bestTune$.C)
    
    if(nlevels(grp) == 2){
        svm.index <- svmrfeFeatureRanking(dat, grp, best.C)
    }else{
        svm.index <- svmrfeFeatureRankingForMulticlass(dat, grp, best.C)
    }
    
    if(model.features){
        warning("SVM currently doesn't have an internal metric or 
                general criteria for optimal number of features.
                \nTop 10% features returned instead")
        top.10 <- round(length(svm.index)/10, 0)
        out <- colnames(dat[,svm.index])[1:top.10]
    }else{
        mod.features <- colnames(dat[,svm.index])
        orig.names <- colnames(dat)
        if(is.null(f)){
            ranks <- c(1:length(svm.index))
            names(ranks) = as.character(mod.features)
            
            out <- ranks[orig.names[orig.names %in% names(ranks)]]
        }else{
            out <- head(mod.features, n = f)
        }
    }
    
    # collect ranked features
    list(features.selected = out)
},

pam = 
{
    if(model.features){
        pam.features <- try(
            data.frame(
                pamr.listgenes(x[[1]], 
                               dat, 
                               threshold = bestTune$.threshold)),
            silent = TRUE)
        if(class(pam.features)[1] == "try-error"){
            tmp <- matrix("", nrow = 1, ncol = 3)
            blank.names <- paste(levels(grp), "Score", sep = "-")
            colnames(tmp) <- c("id", blank.names)
            mod.features <- tmp
        }else{
            mod.features <- as.character(pam.features[,1])
        }
        #}
    }else{
        
        mod.features <- 
            lapply(x, 
                   FUN = function(x){
                       as.character(
                           data.frame(
                               pamr.listgenes(x, dat, threshold = 0)
                           )[,1])
                   })
        
        if(is.null(f)){
            nc <- nrow(dat$x)
            ranks <- c(1:nc)
            names(ranks) <- as.character(mod.features[[1]])
            
            mod.features <- 
                as.data.frame(ranks[dat$geneid[dat$geneid 
                                               %in% names(ranks)]])
            
        }else{
            mod.features <- 
                head(
                    unlist(mod.features, 
                           recursive = FALSE, 
                           use.names = FALSE), 
                    n = f)
        }
    }                    
    
    # collect ranked features
    list(features.selected = mod.features)
},         

glmnet = 
{
    # extract coefficients and remove intercept
    if(nlevels(grp) > 2){
        if(model.features){
            lambda <- bestTune$.lambda
            full.coefs <- coef(x[[1]], s = lambda)
            full.coefs <- 
                lapply(full.coefs, 
                       FUN = function(x){
                           as.matrix(x[2:(ncol(dat)+1),, drop = FALSE])
                       })
            coefs <- 
                lapply(full.coefs, 
                       FUN = function(x) x[x[,1] >= 0,, drop = FALSE])
            
            coef.names <- lapply(coefs, row.names)
            coefs <- unlist(coefs, use.names = TRUE)
            names(coefs) <- unlist(coef.names)
            
            
            var.names <- unique(names(coefs))
            
            for(n in seq(along = unique(names(coefs)))){
                ind <- which(names(coefs) == var.names[n])
                uni <- sum(abs(coefs[ind]))
                names(uni) <- var.names[n]
                coefs <- coefs[-ind]
                coefs <- c(coefs, uni)
            }        
            
            coefs <- coefs[which(coefs > 0)]
            coefs <- coefs[order(-coefs)]
            mod.features <- names(coefs)                            
        }else{
            if(is.null(f)){
                full.coefs <- coef(x[[1]], s = 0)
                full.coefs <- 
                    lapply(full.coefs, 
                           FUN = function(x){
                               as.matrix(x[2:(ncol(dat)+1),, drop = FALSE])
                           })
                coefs <- 
                    lapply(full.coefs, 
                           FUN = function(x) x[x[,1] >= 0,, drop = FALSE])
                
                coef.names <- lapply(coefs, row.names)
                coefs <- unlist(coefs, use.names = TRUE)
                names(coefs) <- unlist(coef.names)
                
                var.names <- unique(names(coefs))                  
                
                for(n in seq(along = unique(names(coefs)))){
                    ind <- which(names(coefs) == var.names[n])
                    uni <- sum(abs(coefs[ind]))
                    names(uni) <- var.names[n]
                    coefs <- coefs[-ind]
                    coefs <- c(coefs, uni)
                }   
                
                coefs <- sort(coefs, decreasing = TRUE)
                ranks <- seq(length(names(coefs)))
                names(ranks) <- names(coefs)
                orig.order <- x[[1]]$xNames
                
                mod.features <- as.data.frame(ranks)[orig.order,, drop = FALSE]
                
                colnames(mod.features) <- "glmnet"
                
            }else{
                # check if a penalized model exists which 
                # includes all features
                index <- lapply(x, FUN = function(x) min(which(x$df >= f)))
                lambda <- 
                    mapply(x, 
                           FUN = function(x, ind) x$lambda[ind], ind = index)
                
                full.coefs <- coef(x[[1]], s = lambda)
                full.coefs <- 
                    lapply(full.coefs, 
                           FUN = function(x){
                               as.matrix(x[2:(ncol(dat)+1),, drop = FALSE])
                           })
                coefs <- 
                    lapply(full.coefs, 
                           FUN = function(x) x[x[,1] != 0,, drop = FALSE])
                
                coef.names <- lapply(coefs, row.names)
                coefs <- unlist(coefs, use.names = TRUE)
                names(coefs) <- unlist(coef.names)
                
                dups <- table(names(coefs))
                dups <- names(dups[dups > 1])
                
                for(n in seq(along = dups)){
                    ind <- which(names(coefs) == dups[n])
                    uni <- sum(abs(coefs[ind]))
                    names(uni) <- dups[n]
                    coefs <- coefs[-ind]
                    coefs <- c(coefs, uni)
                }
                
                mod.features <- 
                    as.data.frame(
                        names(
                            head(
                                abs(coefs)[order(-abs(coefs))], 
                                f)
                        ))
                colnames(mod.features) <- "glmnet"
                
            } # end of f = NULL
        } # end of multinomial feature extraction
    }else{              
        if(model.features){
            lambda <- bestTune$.lambda
            
            if(length(lambda) == 1){
                coefficients <- 
                    mapply(x = x, 
                           FUN = function(x, y){
                               as.matrix(
                                   abs(
                                       coef(x, s = lambda)[2:(ncol(y)+1),
                                                           ,drop=FALSE]
                                   )
                               )
                           },y = list(dat))
                rownames(coefficients) <- colnames(dat)
            }else{
                coefficients <- 
                    mapply(x = x, 
                           FUN = function(x, y, lamb){
                               as.matrix(
                                   abs(
                                       coef(x, s = lamb)[2:(ncol(y)+1),,
                                                         drop=FALSE]
                                   )
                               )
                           }, 
                           y = list(dat), 
                           lamb = list(lambda))
                rownames(coefficients) <- colnames(dat)
            }
            nonzero.coefficients <- 
                apply(coefficients, 2, FUN = function(x) x[x!=0])
            mod.features <- 
                rownames(
                    nonzero.coefficients[order(-nonzero.coefficients),
                                         ,drop = FALSE])
            
        }else{
            if(is.null(f)){
                mod.features <- 
                    mapply(x, 
                           FUN = function(x, y){
                               rank(
                                   -abs(
                                       coef(x, s = 0)[2:(ncol(y)+1),]
                                   )
                               )
                           }, 
                           y = list(dat))
                mod.features <- as.list(as.data.frame(mod.features))
                
                orig.names <- 
                    lapply(mod.features, 
                           FUN = function(x) names(x) = colnames(dat))
                for(i in seq(along = mod.features)){
                    names(mod.features[[i]]) <- orig.names[[i]]
                }
            }else{
                # check if a penalized model exists which includes all features
                index <- 
                    lapply(x, FUN = function(x) min(which(x$df >= f)))
                lambda <- 
                    mapply(x, 
                           FUN = function(x, ind) x$lambda[ind], 
                           ind = index)
                
                if(length(lambda) == 1){
                    coefficients <- 
                        mapply(x, 
                               FUN = function(x, y){
                                   as.matrix(
                                       abs(
                                           coef(x, s = lambda)[2:(ncol(y)+1),
                                                               ,drop=FALSE]
                                       )
                                   )
                               }, 
                               y = list(dat))
                    rownames(coefficients) <- colnames(dat)
                }else{
                    coefficients <- 
                        mapply(x, 
                               FUN = function(x, lamb, dat){
                                   as.matrix(
                                       coef(x, s = lamb)[2:(ncol(dat)+1),
                                                         ,drop=FALSE])
                               }, 
                               lamb = list(lambda), 
                               dat = list(dat))
                    rownames(coefficients) <- colnames(dat)
                }
                
                nonzero.coefficients <- 
                    apply(coefficients, 2, FUN = function(x) x[x!=0])
                mod.features <- 
                    head(
                        rownames(
                            nonzero.coefficients[order(-nonzero.coefficients),
                                                 ,drop = FALSE]), 
                        f)
            }
        }
    }
    
    # collect ranked features
    list(features.selected = mod.features)
}
            )
features
    }
