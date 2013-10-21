# Verify inputs for analysis

# x: matrix or data frame with explanatory variables
# y: vector or factor with group memberships
# na.rm: logical indicating missing values in x

my_verify <- function(x, y, na.rm=na.rm)
  {
    # x matrix or data.frame
    if (is.null(dim(x))) {
      stop("\n Error: 'variables' is not a matrix or dataframe")
    }
    # check lengths of x and y
    if (nrow(x) != length(y))
      stop("\n Error: 'variables' and 'group' have different lengths")
    # no missing values allowed when na.rm=FALSE
    if (!na.rm) {
      if (length(complete.cases(x)) != nrow(x))
        stop("\n Error: no missing values allowed in 'variables'")    
    }
    # y vector or factor
    if (!is.vector(y) && !is.factor(y))
      stop("\n Error: 'group' must be a factor")
    # make sure y is a factor
    if (!is.factor(y)) y = as.factor(y)
    # no missing values in y
    if (any(!is.finite(y)))
      stop("\n Error: no missing values allowed in 'group'")
    # quantitative data
    # make sure is matrix
    if (!is.matrix(x)) x <- as.matrix(x)
    # only numeric values
    if (!is.numeric(x))
      stop("\n Error: 'variables' must contain only numeric values")
    
    # verified inputs
    # fixed colnames generation
    # was rep("X",1:ncol(x),sep='')
    if (is.null(colnames(x))) colnames(x) = paste(rep("X",ncol(x)),seq_len(ncol(x)), sep='') 
    if (is.null(rownames(x))) rownames(x) = 1:nrow(x)
    list(X=x, Y=y)
}
