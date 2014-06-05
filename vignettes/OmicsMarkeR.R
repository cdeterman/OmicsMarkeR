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
      )),
  D = 10
  ## how many discriminating variables
  )$discr.mat


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
                        f = 10, 
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
results$RPT

###################################################
### code chunk number 6: Alternate Tuning Resolutions
###################################################
grid <- denovo.grid(vars, "plsda", 15)
grid <- append(grid, denovo.grid(vars, "svm", 3)) 

# list of two different methods resolution

results <- fs.stability(vars, 
                        groups, 
                        grid = grid,
                        method = c("plsda", "svm"), 
                        f = 10, 
                        k = 3, 
                        k.folds = 10, 
                        verbose = FALSE)

###################################################
### code chunk number 7: Manual Tuning Resolutions
###################################################
grid <- list(expand.grid(.ncomp = seq(3,8)))
grid <- append(grid, list(expand.grid(.C = c(.50, 1))))
names(grid) <- c("plsda", "svm")