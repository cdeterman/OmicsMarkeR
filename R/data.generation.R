###################################
### Simulated Metabolomics Data ###
###################################

plus_minus <- function(beta){
  inv <-rbinom(beta,1,0.5)
  inv <-ifelse(inv==0,-1,1)
  return(inv)
}

ind_corr <- function(dat, group_size, start, end){
  pm <- plus_minus(group_size-1)
  
  # take first column of block
  f <- dat[,start]
  
  #copy to each subsequent column
  dat[,(start+1):end] <- f
  dat[,(start+1):end]

  #change sign
  finish <- t(t(dat[,(start+1):end])*pm)  
  return(finish)
}

#' @title Noise Matrix Generator
#' @description Provides a matrix to perturb randomly generated data to facilitate a more
#' realistic dataset.
#' @param matrix A matrix of simulated data with dimensions comparable to 'real' datasets
#' @param k Correlation Perturbation - The higher k, the more the data is perturbed.
#' @return Returns a matrix of the same dimensions as \code{matrix} that can add to perturb the 
#' original simulated data.
#' @author Charles E. Determan Jr.
#' @export

noise.matrix <- function(matrix, k){
  nvar <- ncol(matrix)
  nsamp <- nrow(matrix)
  
  # for each variable j a value wj is obtained from a random generated uniform distribution from 0 to k
  wj <- runif(nvar, max=k, min=0)
  
  # wij - vector of uniform random numbers -Wj < wij < Wj
  # add wij to each element of matrix B
  noise <- matrix(0, nrow=nsamp, ncol=nvar)
  for (i in 1:ncol(noise)){
    wij <- runif(nsamp, min=-wj[i], max=wj[i])
    noise[,i] <- wij}
  
  return(noise)
}

#' @title Random Multivariate Data Generator
#' @description Generates a matrix of dimensions \code{nvar} by \code{nsamp} consisting of 
#' random numbers generated from a normal distriubtion.  This normal distribution is then
#' perturbed to more accurately reflect experimentally acquired multivariate data.
#' @param nvar Number of features (i.e. variables)
#' @param nsamp Number of samples
#' @param st.dev The variation (i.e. standard deviation) that is typical in datasets 
#' of interest to the user.  Default \code{spread = 1}
#' @param perturb The amount of perturbation to the normal distribution.  Default \code{perturb = 0.2}
#' @return Matrix of dimension \code{nvar} by \code{nsamp}
#' @author Charles E. Determan Jr.
#' @seealso \code{\link{create.corr.matrix}}, \code{\link{create.discr.matrix}}
#' @export

create.random.matrix <-
  function(nvar,
           nsamp,
           st.dev = 1,
           perturb = 0.2
           )
    {
    ## Random numbers matrix
    R <- replicate(nvar, rnorm(nsamp, mean=0, sd=st.dev))
    
    ### Perturb Normal Distribution
    # Random uniform numbers matrix
    N <- replicate(nvar, runif(nsamp, max=perturb, min=-perturb))
    
    # Null matrix w/o any correlations
    U <- R + N
    U
  }

#' @title Correlated Multivariate Data Generator
#' @description Generates a matrix of dimensions \code{dim(U)} with induced correlations.  Blocks
#' of variables are randomly assigned and correlations are induced.  A noise matrix is applied
#' to the final matrix to perturb 'perfect' correlations.
#' @param U Numeric matrix
#' @param k Correlation Perturbation - The higher k, the more the data is perturbed. Default \code{k = 4}
#' @param min.block.size minimum number of variables to correlate Default \code{min.block.size = 2}
#' @param max.block.size maximum number of variables to correlate Default \code{max.block.size = 5}
#' @return Matrix of dimension \code{dim(U)} with correlations induced between variables
#' @note Output does not contain classes, may provide externally as classes are irrelevant in this function.
#' @author Charles E. Determan Jr.
#' @seealso \code{\link{create.random.matrix}}, \code{\link{create.discr.matrix}}
#' @export

create.corr.matrix <-
  function(U,
           k=4,
           min.block.size = 2,
           max.block.size = 5
           )
    {
    if(any(!sapply(U, is.numeric))){
      stop("Error: matrix components must all be numeric.  Check to make sure no factors are in matrix.")
    }
    
    nvar <- ncol(U)
    nsamp <- nrow(U)
        
    # generate block size (groups of variables to correlate)
    # discrete uniform distribution, 2>beta>5 - original (more appropriate for MS?)
    # set beta -> 1>x>5 to reflect that some variables likely independent with less dimensions in NMR
    beta <- sample(min.block.size:max.block.size, size=nvar-1, replace=T)    
    
    # loop through each block
    for (i in 1:length(beta)){ 
      if (i == 1){
        start <- 1
        end <- start + beta[i] - 1
        
        U[,(start+1):end] <- ind_corr(U, beta[i], start, end)
      }else{
        start <- end + 1
        
        if(start > nvar){
          start <- nvar
          end <- nvar
          cat('solo last variable')
          break
        }
        
        end <- end + beta[i]
        
        if(end > nvar){
          end <- nvar
          is_group <- end - start
          if (is_group > 1){
            U[,(start+1):end] <- ind_corr(U, is_group, start, end)
          }
          break
        }
        
        if (beta[i] != 1){
          U[,(start+1):end] <- ind_corr(U, beta[i], start, end)  
        }  
        B <- U
      }
    }
    
    
    ## Create a noise matrix to perturb distributions
    W <- noise.matrix(B, k)
    B <- B + W
    
    # randomize order of columns to make blocks no longer obvious
    rand <- sample(nsamp)
    V <- structure(B[rand,], class = "matrix.corr")
    V
  }

#' @title Discriminatory Multivariate Data Generator
#' @description Generates a matrix of dimensions \code{dim(U)} with induced correlations.  \code{D}
#' variables are randomly selected as discriminatory.  If \code{num.groups = 2} then discrimination
#' is induced by adding and subtracting values derived from the level of of discrimination, \code{l},
#' for the classes respectively.  Multi-class datasets have a few further levels of randomization.  For
#' each variable, a random number of the groups are selected as discriminating while the remaining groups
#' are not altered.  For each discriminatory group, a unique change is provided by randomly assigning
#' addition or subtraction of the discrimination factor.  For example, if 3 groups are selected and two 
#' groups are assigned as addition and the third subtraction, the second addition is multiplied by its
#' number of replicates.  E.g. (1,1,-1) -> (1,2,-1).  These values are randomized and then multiplied by
#' the respective discrimination factor.  The resulting values are then added/subtracted from the 
#' respective groups.  A noise matrix is applied to the final matrix to perturb 'perfect' 
#' discrimination.
#' @param V Numeric matrix of class "matrix.corr"
#' @param D Number of discriminatory variables induced. Default \code{D = 20}
#' @param l Level of discrimination, higher = greater separation.  Default \code{l = 1.5}
#' @param num.groups Number of groups in the dataset
#' @param k Correlation Perturbation - The higher k, the more the data is perturbed. Default \code{k = 4}
#' @return Matrix of dimension \code{dim(V)+1} with discriminatory variables induced and the 
#' .classes added to the end of the matrix.
#' @author Charles E. Determan Jr.
#' @import data.table
#' @export

create.discr.matrix <-
  function(V,
           D = 20,
           l = 1.5,
           num.groups = 2,
           k = 4
           )
  {
    if(class(V) != "matrix.corr"){
      warning("Matrix provided was not produced from create.corr.matrix. Assumptions are unknown.  Matrix may assume complete independence and/or normality.  Proceed with caution")
    }
    
    nc <- ncol(V)
    nr <- nrow(V)
    groups <- LETTERS[seq(num.groups)]
    classes <- unlist(mapply(groups, FUN = function(x,y) rep(x, nrow(y)), y = Z.list), use.names=FALSE)    
    
    # randomly select which variables to be discriminatory
    d <- sample(nc, size = D, replace=F)
    Z <- V[,d]
        
    # get range of discriminatory ability (as is typical in real data)
    di <- runif(D, min=-l, max=l)
    
    if(num.groups == 2){
      # Binary class induced discrimination
      I <- nr/2
      # add di to first group, subtract from second to induce discrimination
      for(i in seq(D)){
        Z[1:I,i] <- Z[1:I,i]+di[i]
        Z[(1+I):nr,i] <- Z[(1+I):nr,i]-di[i] 
      }
    }else{
      # Multi-class induced discrimination
      # split Z into n groups
      Z.list <- split(as.data.frame(Z), rep(1:num.groups, each = nrow(Z)/num.groups))
      Z.rows <- lapply(Z.list, rownames)
      
      # Create matrix of instructions for discriminating variables
      mat <- matrix(0, nrow=D, ncol=num.groups)
      rownames(mat) <- paste("Var", seq(D), sep=".")
      for(d in seq(D)){
        # for each variable, randomize number of groups discriminating
        num <- sample.int(num.groups, 1)
        
        # determine how values will be added/subtracted
        pm <- plus_minus(num)
        
        # remainder of groups retain same value
        num.zero <- num.groups - num
        pm.complete <- c(pm, rep(0, num.zero))
        
        # randomize how value added
        mat[d,] <- sample(pm.complete)
      }
      
      # make each change unique to facilitate each group is discriminated while leaving 0's alone
      mat.uniq <- t(apply(mat, 1, FUN = function(m) unsplit(lapply(split(m, factor(m)), FUN = function(x) x*seq(length(x))), factor(m))))
      # multiple direction changes by di value
      mat.discr <- mat.uniq*di
      
      for(l in seq(along = Z.list)){
        # randomize order of columns so different groups are changed differently
        mat.rand <- mat.discr[,sample(ncol(mat.discr))]
        Z.dat <- Z.list[[l]]
        #vector for each group
        di.sub <- mat.rand[,l]
        Z.list[[l]] <- sweep(Z.dat, MARGIN = 2, di.sub, '+')
        rownames(Z.list[[l]]) <- Z.rows[[l]]
      }
      
    }
    #require(data.table)
    Z <- rbindlist(Z.list)
    
    # add newly discriminated variables back to matrix V to make matrix S
    S <- V
    S[,d] <- Z
    
    ## Create a noise matrix to perturb distributions again
    W <- noise.matrix(S, k)
    Y <- S + W
    Y$.classes <- classes
    Y
  }




# histograms of correlation coefficients

# function to determine max correlation
# must be the second highest because always 1 for correlation with self
# therefore remove max in each column and return new max
#f4 <- function (x)
#  {
#  apply(x, 2, function(column) max(column[-which.max(column)]))
#}

#hist(f4(abs(cor(U))))
#hist(f4(abs(cor(B))))
#hist(f4(abs(cor(Y))))

