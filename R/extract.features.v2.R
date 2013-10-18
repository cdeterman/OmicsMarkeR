#' @title Feature Extraction
#' @description Extracts features from models that have been previously fit.
#' @param x Previously fitted model
#' @param dat Numeric variable data used for fitted models (In appropriate format)
#' @param method String indicating the INDIVIDUAL model being extracted from
#' @param model.features Logical argument dictating if features selected determined by models instead of
#' user determined number of features.
#' @param bestTune If \code{model.features = TRUE}, must provide the parameter at which to extract
#' features from the model.
#' @param comp.catch An internal check for plsda models.  If the optimal model contains only 1 component,
#' the ncomp paramter must be set to 2 for the model.  However, features are still extracted only from the first component.
#' @return Returns list of the features selected from the fitted model.
#' @export

extract.features <-
  function(x,                       # fitted model
           dat = NULL,              # training variable data or pam format
           grp = NULL,              # training class data
           method,                  # which model being fit
           model.features = FALSE,  # state if model defined variable selection performed
           bestTune = NULL,         # tune metric for svm and pam if model.features = TRUE
           f,                       # number of features to subset if defined
           comp.catch = NULL        # a check for plsda models that have only 1 component
  )
{    
    features <- switch(class(x[[1]]),
                       plsda =
                         {
                           if(!is.null(comp.catch)){
                             l <- 1
                           }else{
                             l <- ncol(x[[1]]$VIP)  
                           }
                           
                           if(model.features){
                             # generally a VIP >= 1 is significant in the literature
                             sub <- lapply(x, FUN = function(x) x$VIP[which(x$VIP[,l] >= 1),])
                             mod.features <- lapply(sub, FUN = function(x) names(sort(x[,l], decreasing=T)))
                           }else{
                             if(is.null(f)){
                               # plsda VIP = higher is more important so must reverse rank function with '-'
                               mod.features <- lapply(x, FUN = function(x) rank(-x$VIP[,l]))
                             }else{
                               mod.features <- lapply(x, FUN = function(x) names(sort(x$VIP[,l], decreasing = T)))
                               mod.features <- lapply(mod.features, head, n = f)
                             }
                           }
                           
                           ## Make a matrix from the list, with shorter vectors filled out with ""
                           n <- max(sapply(mod.features, FUN = function(x) length(x)))
                           ll <- lapply(mod.features, function(x) {
                             c(as.character(x), rep("", times = n - length(x)))
                             #c(as.character(X), rep("", times = testn - length(X)))
                           })
                           out <- as.data.frame(do.call(cbind, ll)) 
                           
                           # collect ranked features
                           list(features.selected = out)
                         },
                       
                       gbm =  
                         {
                           if(model.features){
                             # generally a VIP >= 1 is significant in the literature
                             vips <- lapply(x, FUN = function(x) relative.influence(x, n.trees = bestTune$.n.trees))
                             sub <- lapply(vips, FUN = function(x) x[which(x >= 1)])
                             mod.features <- lapply(sub, FUN = function(x) names(sort(x, decreasing=T)))
                           }else{
                             if(is.null(f)){
                               # plsda VIP = higher is more important so must reverse rank function with '-'
                               mod.features <- lapply(x, FUN = function(x) rank(-relative.influence(x, n.trees = x$n.trees)))
                               }else{
                                 mod.features <- lapply(x, FUN = function(x) names(relative.influence(x, n.trees = x$n.trees, sort=T)))
                                 mod.features <- lapply(mod.features, head, n = f)
                             }
                           }
                           
                           ## Make a matrix from the list, with shorter vectors filled out with ""
                           n <- max(sapply(mod.features, FUN = function(x) length(x)))
                           ll <- lapply(mod.features, function(x) {
                             c(as.character(x), rep("", times = n - length(x)))
                             #c(as.character(X), rep("", times = testn - length(X)))
                           })
                           out <- as.data.frame(do.call(cbind, ll)) 
                           
                           # collect ranked features
                           list(features.selected = out)
                         },
                       
                       randomForest =
                         {
                           # Mean Decrease in Accuracy metric (type=1)
                           # Gini Index (type=2)
                           
                           if(model.features){
                             # generally a VIP >= 1 is significant in the literature
                             vips <- lapply(x, FUN = function(x) importance(x, type = 1))
                             sub <- lapply(vips, FUN = function(x) x[which(x >= 1),,drop = F])
                             mod.features <- lapply(sub, FUN = function(x) names(x[sort(x, decreasing = T, index.return= T)$ix,]))
                           }else{
                             if(is.null(f)){
                               # random forest VIP = higher is more important so must reverse rank function with '-'
                               mod.features <- lapply(x, FUN = function(x) rank(-importance(x, type = 1)))
                             }else{
                               vips <- lapply(x, FUN = function(x) importance(x, type = 1))
                               mod.features <- lapply(vips, FUN = function(x) names(x[sort(x, decreasing = T, index.return=T)$ix,]))
                               mod.features <- lapply(mod.features, head, n = f)
                             }
                           }
                           
                           ## Make a matrix from the list, with shorter vectors filled out with ""
                           n <- max(sapply(mod.features, FUN = function(x) length(x)))
                           ll <- lapply(mod.features, function(x) {
                             c(as.character(x), rep("", times = n - length(x)))
                             #c(as.character(X), rep("", times = testn - length(X)))
                           })
                           out <- as.data.frame(do.call(cbind, ll)) 
                           
                           # collect ranked features
                           list(features.selected = out)
                         },                  
                       
                       svm =   
                         {
                           best.C <- vector("list", length(dat))
                           for(i in seq(along = best.C)){
                             best.C[[i]] <- bestTune$.C
                           }
                           
                           if(nlevels(grp) == 2){
                             best.C <- vector("list", length(dat))
                             for(i in seq(along = best.C)){
                               best.C[[i]] <- bestTune$.C
                             }
                             svm.index <- mapply(dat, FUN = function(x,y,z) svmrfeFeatureRanking(x, y, z), y = grp, z = best.C)
                           }else{
                             svm.index <- mapply(dat, FUN = function(x,y,z) svmrfeFeatureRankingMulticlass(x, y, z), y = grp, z = best.C)
                           }
                           
                           if(model.features){
                             warning("SVM currently doesn't have an internal metric or general criteria for optimal number of features.\nTop 10% features returned instead")
                             top.10 <- round(nrow(svm.index)/10, 0)
                             out <- apply(svm.index, 2, FUN = function(x) colnames(dat[[1]][,x])[1:top.10])
                           }else{
                             mod.features <- apply(svm.index, 2, FUN = function(x) colnames(dat[[1]][,x]))
                             orig.names <- colnames(dat[[1]])
                             if(is.null(f)){
                               ranks <- rep(list(ranks = 1:nrow(svm.index)), 5)
                               for(i in seq(along = ranks)){
                                 names(ranks[[i]]) <- as.character(mod.features[,i])
                               }
                               out <- lapply(ranks, FUN = function(x) x[orig.names[orig.names %in% names(x)]])
                               }else{
                                 out <- head(mod.features, n = f)
                               }
                           }
                           
                           # collect ranked features
                           list(features.selected = out)
                         },
                       
                       pamrtrained = 
                         {               
                           #best.threshold <- vector("list", length(dat))
                           #for(i in seq(along = best.C)){
                          #   best.threshold[[i]] <- bestTune$.threshold
                          # }
                          
                           
                           #data.frame(pamr.listgenes(x[[2]], dat[[2]], threshold = best.threshold[[2]]))
                           #x <- finalModel
                           
                           if(model.features){
                             mod.features <- vector("list", length(x))
                             for(i in seq(along = x)){
                               pam.features <- try(mod.features[[i]] <- data.frame(pamr.listgenes(x[[i]], 
                                                                                                  dat[[i]], 
                                                                                                  threshold = bestTune$.threshold)),
                                                   silent = TRUE)
                               if(class(pam.features)[1] == "try-error"){
                                 tmp <- matrix("", nrow = 1, ncol = 3)
                                 colnames(tmp) <- c("id", "A-Score", "B-Score")
                                 mod.features[[i]] <- tmp
                               }
                             }
                             }else{
                               for(i in seq(along = dat)){
                                 mod.features <- lapply(x, FUN = function(x) data.frame(pamr.listgenes(x, dat[[i]], threshold = 0)))
                               }
                             
                             if(is.null(f)){
                               
                               ranks <- rep(list(ranks = 1:nc), length(x))
                               for(i in seq(along = ranks)){
                                 names(ranks[[i]]) <- as.character(mod.features[[i]][,1])
                               }
                               
                               for(i in seq(along = ranks)){
                                 mod.features <- lapply(ranks, FUN = function(x) as.data.frame(x[dat[[i]]$geneid[dat[[i]]$geneid %in% names(x)]]))
                               }

                             }else{
                               mod.features <- lapply(mod.features, head, n = f)
                             }
                           }
                           
                           ## Make a matrix from the list, with shorter vectors filled out with ""
                           n <- max(sapply(mod.features, FUN = function(x) length(x[,1])))
                           ll <- lapply(mod.features, function(x) {
                             c(as.character(x[,1]), rep("", times = n - length(x[,1])))
                             #c(as.character(X), rep("", times = testn - length(X)))
                           })
                           out <- as.data.frame(do.call(cbind, ll))                           
                           
                           # collect ranked features
                           #list(features.selected = rfs)
                           list(features.selected = out)
                         },         
                       
                       glmnet = 
                        {
                          
                          x <- finalModel
                          
                          # extract coefficients and remove intercept
                          if(nlevels(dat$.classes) > 2){
                            if(model.features){
                              
                              #coefficients <- as.matrix(abs(coef(x, s = bestTune$.lambdaOpt)[[1]][2:(nc+1),,drop=FALSE]))
                              lambda <- unlist(lapply(x, FUN = function(x) x$lambdaOpt))
                              
                              if(length(lambda) == 1){
                                mapply(x = x, FUN = function(x, y) as.matrix(abs(coef(x, s = lambda)[[1]][2:(ncol(y)+1),,drop=FALSE])), y = dat)
                              }
                              #str(x[[1]])
                              coefficients <- lapply(x, FUN = function(x) as.matrix(abs(coef(x, s = x$lambdaOpt)[[1]][2:(ncol(dat)),,drop=FALSE])))
                              mod.features <- lapply(coefficients, FUN = function(x) rownames(x[order(-x),,drop=F]))
                              
                            }else{
                              if(is.null(f)){
                                #rfs[,iter] <- rank(-abs(coef(x, s = 0)[[1]][2:(nc+1),]))
                                mod.features <- lapply(x, FUN = function(x) rank(-abs(coef(x, s = 0)[[1]][2:ncol(dat),])))
                                
                                }else{
                                  # check if a penalized model exists which includes all features
                                  #index <- min(which(x$df >= f))
                                  #lambda <- x$lambda[index]
                                  index <- lapply(x, FUN = function(x) min(which(x$df >= f)))
                                  lambda <- mapply(x, FUN = function(x, ind) x$lambda[ind], ind = index)
                                  
                                  #rfs[,iter] <- head(rownames(coefficients[order(-coefficients),, drop=F]), f)
                                  #coefficients <- mapply(x, FUN = function(x, lamb, dat) as.matrix(abs(coef(x, s = lamb)[[1]][2:(ncol(dat)),,drop=FALSE])), lamb = lambda, dat = dat)
                                  mapply(x, FUN = function(x, lamb, dat) as.matrix(coef(x, s = lamb)[[1]][2:(ncol(dat)),,drop=FALSE]), lamb = lambda, dat = dat)
                                  
                                  #coefficients <- as.matrix(abs(coef(x, s = bestTune$.lambdaOpt)[[1]][2:(nc+1),,drop=FALSE]))
                                  
                                  as.matrix(abs(coef(x[[1]], s = as.numeric(lambda[1]))[[1]][2:(ncol(dat[[1]])),,drop=FALSE]))
                                  
                                  #ncol(dat[[1]])
                                  mod.features <- lapply(coefficients, FUN = function(x) head(rownames(x[order(-x),,drop=F]),f))
                              }
                              
                            }
                            # end of multinomial feature extraction
                          }else{
                              #coefficients <- as.matrix(abs(coef(x, s = lambda)[2:(nc+1),,drop=FALSE]))
                              
                              if(model.features){
                                #coefficients <- as.matrix(abs(coef(x, s = bestTune$.lambdaOpt)[[1]][2:(nc+1),,drop=FALSE]))
                                lambda <- unlist(lapply(x, FUN = function(x) x$lambdaOpt))
                                
                                if(length(lambda) == 1){
                                  coefficients <- mapply(x = x, FUN = function(x, y) as.matrix(abs(coef(x, s = lambda)[2:(ncol(y)+1),,drop=FALSE])), y = dat)
                                  rownames(coefficients) <- colnames(dat[[1]])
                                }else{
                                  coefficients <- mapply(x = x, FUN = function(x, y, lamb) as.matrix(abs(coef(x, s = lamb)[2:(ncol(y)+1),,drop=FALSE])), y = dat, lamb = lambda)
                                  rownames(coefficients) <- colnames(dat[[1]])
                                }
                                #str(x[[1]])
                                #coefficients <- lapply(x, FUN = function(x) as.matrix(abs(coef(x, s = x$lambdaOpt)[2:(ncol(dat)),,drop=FALSE])))
                                nonzero.coefficients <- apply(coefficients, 2, FUN = function(x) x[x!=0])
                                #tmp
                                mod.features <- lapply(nonzero.coefficients, FUN = function(x) names(x[order(-x)]))
                                #mod.features <- do.call("cbind", tmp)
                                
                                }else{
                                  if(is.null(f)){
                                    #rfs[,iter] <- rank(-abs(coef(x, s = 0)[[1]][2:(nc+1),]))
                                    #mod.features <- lapply(x, FUN = function(x) rank(-abs(coef(x, s = 0)[2:ncol(dat),])))
                                    mod.features <- mapply(x, FUN = function(x, y) rank(-abs(coef(x, s = 0)[2:(ncol(y)+1),])), y = dat)
                                    mod.features <- as.list(as.data.frame(mod.features))
                                    
                                    orig.names <- lapply(mod.features, FUN = function(x) names(x) = colnames(dat[[1]]))
                                    for(i in seq(along = mod.features)){
                                      names(mod.features[[i]]) <- orig.names[[i]]
                                    }
                                                                        
                                    
                                    
                                  }else{
                                    # check if a penalized model exists which includes all features
                                    #index <- min(which(x$df >= f))
                                    #lambda <- x$lambda[index]
                                    index <- lapply(x, FUN = function(x) min(which(x$df >= f)))
                                    lambda <- mapply(x, FUN = function(x, ind) x$lambda[ind], ind = index)
                                    
                                    #rfs[,iter] <- head(rownames(coefficients[order(-coefficients),, drop=F]), f)
                                    #coefficients <- mapply(x, FUN = function(x, lamb, dat) as.matrix(abs(coef(x, s = lamb)[[1]][2:(ncol(dat)),,drop=FALSE])), lamb = lambda, dat = dat)
                                    
                                    if(length(lambda) == 1){
                                      coefficients <- mapply(x, FUN = function(x, y) as.matrix(abs(coef(x, s = lambda)[2:(ncol(y)+1),,drop=FALSE])), y = dat)
                                      rownames(coefficients) <- colnames(dat[[1]])
                                    }else{
                                      coefficients <- mapply(x, FUN = function(x, lamb, dat) as.matrix(coef(x, s = lamb)[2:(ncol(dat)+1),,drop=FALSE]), lamb = lambda, dat = dat)
                                      rownames(coefficients) <- colnames(dat[[1]])
                                    }
                                    
                                    
                                    
                                    nonzero.coefficients <- apply(coefficients, 2, FUN = function(x) x[x!=0])
                                    #tmp
                                    mod.features <- lapply(nonzero.coefficients, FUN = function(x) head(names(x[order(-x)]), f))
                                  }
                                  
                                }
                          }
                          
                          
                          ## Make a matrix from the list, with shorter vectors filled out with ""
                          n <- max(sapply(mod.features, FUN = function(x) length(x)))
                          ll <- lapply(mod.features, function(x) {
                            c(as.character(x), rep("", times = n - length(x)))
                            #c(as.character(X), rep("", times = testn - length(X)))
                          })
                          out <- as.data.frame(do.call(cbind, ll))     
                          
                          
                          # collect ranked features
                          #list(features.selected = rfs)
                          list(features.selected = out)
                        }
    )
    features
  }
