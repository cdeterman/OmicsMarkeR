## ----style, eval=TRUE, echo=FALSE, results="asis"---------------------------------------
BiocStyle::latex()

## ----install, eval = FALSE--------------------------------------------------------------
#  source("http://bioconductor.org/biocLite.R")
#  biocLite("OmicsMarkeR")

## ----datagen----------------------------------------------------------------------------
library("OmicsMarkeR")
set.seed(123)
dat.discr <- create.discr.matrix(
    create.corr.matrix(
        create.random.matrix(nvar = 50, 
                             nsamp = 100, 
                             st.dev = 1, 
                             perturb = 0.2)),
    D = 10
)

## ----fs.stability-----------------------------------------------------------------------
vars <- dat.discr$discr.mat
groups <- dat.discr$classes
fits <- fs.stability(vars, 
                    groups, 
                    method = c("plsda", "rf"), 
                    f = 10, 
                    k = 3, 
                    k.folds = 10, 
                    verbose = 'none')

## ----performance------------------------------------------------------------------------
performance.metrics(fits)
fits$RPT

## ----feature.table----------------------------------------------------------------------
feature.table(fits, "plsda")

## ----predictClasses, eval=FALSE---------------------------------------------------------
#  
#  # create some 'new' data
#  newdata <- create.discr.matrix(
#      create.corr.matrix(
#          create.random.matrix(nvar = 50,
#                               nsamp = 100,
#                               st.dev = 1,
#                               perturb = 0.2)),
#      D = 10
#  )$discr.mat
#  
#  # original data combined to a data.frame
#  orig.df <- data.frame(vars, groups)
#  
#  # see what the PLSDA predicts for the new data
#  # NOTE, newdata does not require a .classes column
#  predictNewClasses(fits, "plsda", orig.df, newdata)

## ----ensemble, eval=FALSE---------------------------------------------------------------
#  fits <- fs.ensembl.stability(vars,
#                              groups,
#                              method = c("plsda", "rf"),
#                              f = 10,
#                              k = 3,
#                              k.folds = 10,
#                              verbose = 'none')

## ----aggregation------------------------------------------------------------------------
# test data
ranks <- replicate(5, sample(seq(50), 50))
row.names(ranks) <- paste0("V", seq(50))

head(aggregation(ranks, "CLA"))

## ----grid, eval=FALSE-------------------------------------------------------------------
#  # requires data.frame of variables and classes
#  plsda <- denovo.grid(orig.df, "plsda", 3)
#  rf <- denovo.grid(orig.df, "rf", 5)
#  
#  # create grid list
#  # Make sure to assign appropriate model names
#  grid <- list(plsda=plsda, rf=rf)
#  
#  # pass to fs.stability or fs.ensemble.stability
#  fits <- fs.stability(vars,
#                      groups,
#                      method = c("plsda", "rf"),
#                      f = 10,
#                      k = 3,
#                      k.folds = 10,
#                      verbose = 'none',
#                      grid = grid)
#  

## ----metabs-----------------------------------------------------------------------------
metabs <- paste("Metabolite", seq(20), sep="_")

## ----samples----------------------------------------------------------------------------
set.seed(13)
run1 <- sample(metabs, 10)
run2 <- sample(metabs, 10)

## ----jaccard----------------------------------------------------------------------------
jaccard(run1, run2)

## ----kuncheva---------------------------------------------------------------------------
# In this case, 20 original variables
kuncheva(run1, run2, 20)

## ----repeat.metabs----------------------------------------------------------------------
set.seed(21)
# matrix of Metabolites identified (e.g. 5 trials)
features <- replicate(5, sample(metabs, 10))

## ----pairwise.stability-----------------------------------------------------------------
pairwise.stability(features, "sorensen")

## ----model.stability--------------------------------------------------------------------
set.seed(999)
plsda <- 
    replicate(5, paste("Metabolite", sample(metabs, 10), sep="_"))
rf <-
    replicate(5, paste("Metabolite", sample(metabs, 10), sep="_"))

features <- list(plsda=plsda, rf=rf)

# nc may be omitted unless using kuncheva
pairwise.model.stability(features, "kuncheva", nc=20)

## ----permutations, eval=FALSE-----------------------------------------------------------
#  # permuate class
#  perm.class(fits, vars, groups, "rf", k.folds=5,
#             metric="Accuracy", nperm=10)
#  
#  
#  # permute variables/features
#  perm.features(fits, vars, groups, "rf",
#          sig.level = .05, nperm = 10)

## ----doMC, eval=FALSE-------------------------------------------------------------------
#  library(doMC)
#  
#  n <- detectCores()
#  registerDoMC(n)

## ----SNOW, eval=FALSE-------------------------------------------------------------------
#  library(parallel)
#  library(doSNOW)
#  
#  # get number of cores
#  n <- detectCores()
#  
#  # make clusters
#  cl <- makeCluster(n)
#  
#  # register backend
#  registerDoSNOW(cl)

## ----sessionInfo------------------------------------------------------------------------
sessionInfo()

