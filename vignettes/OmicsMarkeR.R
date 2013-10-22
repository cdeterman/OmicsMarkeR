### R code from vignette source 'OmicsMarkeR.Rnw'

###################################################
### code chunk number 1: install (eval = FALSE)
###################################################
## install.packages("OmicsMarkeR", dependencies = c("Imports"))

###################################################
### code chunk number 2: install devel version (eval = FALSE)
###################################################
## library(devtools)
## install_github("OmicsMarkeR", username = "cdeterman")

###################################################
### code chunk number 3: 
###################################################
library("OmicsMarkeR")
set.seed(123)
dat.discr <- create.discr.matrix(
  create.corr.matrix(
    create.random.matrix(
      nvar = 50, 
      ## the number of variables
      nsamp = 100, 
      ## the number of samples
      st.dev = 1, 
      ## the estimated standard deviation
      perturb = 0.2
      ## how noisy the data is
      )))


###################################################
### code chunk number 4: Classification and Feature Selection
###################################################
## isolate variables
vars <- dat.discr[,1:(ncol(dat.discr)-1)]
## isolate groups
groups <- dat.discr[,ncol(dat.discr)]
results <- fs.stability(vars, 
                        groups, 
                        method = c("plsda", "rf"), 
                        f = 20, 
                        ## number of top features returned
                        k = 3, 
                        ## number of bootstrap iterations
                        k.folds = 10, 
                        ## number of k-folds for cross-validation
                        verbose = FALSE
                        ## print output progress
                        )


###################################################
### code chunk number 5: Performance Metrics
###################################################
performance.metrics(results)

