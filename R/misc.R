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


# slight modification of the 'caret' function 'byComplexity' 
# to suit this program
byComplexity2 <- function(x, model)
{
    # must be dataframe to access components easily with '$'
    if(!is.data.frame(x)) x <- as.data.frame(x)
    
    switch(tolower(model),
           gbm =
{
    # As Max Kuhn (caret creator) has stated:
    # This is a toss-up, but the # trees probably adds
    # complexity faster than number of splits
    x[order(x$n.trees, x$interaction.depth, x$shrinkage),] 
},
rf =, plsda =, pam = 
{
    x[order(x[,1]),]
},
svm =
{
    x[order(x$C),]
},
glmnet = x[order(-x$lambda, x$alpha),],
stop("no sorting routine for this model")
    )
}



# These functions or duplicates from the 'caret' package
# They have been copied as they are low-level functions that are not exported
# However, they prove useful for this program as it has a 
# very similar structure
MeanSD <- function(x, exclude = NULL)
{
    if(!is.null(exclude)) x <- x[, !(colnames(x) %in% exclude), drop = FALSE]
    out <- c(colMeans(x, na.rm = TRUE), sapply(x, sd, na.rm = TRUE))
    names(out)[-(1:ncol(x))] <- 
        paste(
            names(out)[-(1:ncol(x))],
            "SD", 
            sep = ""
            )
    out
}

expandParameters <- function(fixed, seq)
{
    if(is.null(seq)) return(fixed)
    
    isSeq <- names(fixed) %in% names(seq)
    out <- fixed
    for(i in 1:nrow(seq))
    {
        tmp <- fixed
        tmp[,isSeq] <- seq[i,]
        out <- rbind(out, tmp)
    }
    out
}

flatTable <- function(pred, obs)
{
    cells <- as.vector(table(pred, obs))
    if(length(cells) == 0) cells <- rep(NA, length(levels(obs))^2)
    names(cells) <- paste(".cell", seq(along= cells), sep = "")
    cells
}

# this is a replication of caret:::progress as it isn't exported
# and the function provides the desired functionality
progress <- function (x, names, iter, start = TRUE) 
{
    text <- paste(ifelse(start, "+ ", "- "), names[iter], ": ", 
                  paste(colnames(x), x, sep = "=", collapse = ", "), sep = "")
    cat(text, "\n")
}

# Options for selecting optimal model by caret
# currently incorporating caret:::best
# Currently don't support oneSE or tolerance
# Dependent upon user needs


cbind.fill <- function(...){
    nm <- list(...) 
    nm <- lapply(nm, as.matrix)
    n <- max(sapply(nm, nrow)) 
    do.call(cbind, lapply(nm, function (x) 
        rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}
