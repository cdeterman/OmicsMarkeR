### R code from vignette source 'OmicsMarkeR.Rnw'

###################################################
### code chunk number 1: install (eval = FALSE)
###################################################
## # currently in development
## # install.packages("OmicsMarkeR", dependencies = c("Imports"))


###################################################
### code chunk number 2: install (eval = FALSE)
###################################################
## library(devtools)
## install_github("OmicsMarkeR", username = "cdeterman")


###################################################
### code chunk number 3: OmicsMarkeR.Rnw:119-127
###################################################
library("OmicsMarkeR")
set.seed(123)
dat.discr <- create.discr.matrix(
  create.corr.matrix(
    create.random.matrix(nvar = 50, 
                         nsamp = 100, 
                         st.dev = 1, 
                         perturb = 0.2)))


###################################################
### code chunk number 4: OmicsMarkeR.Rnw:132-141
###################################################
vars <- dat.discr[,1:(ncol(dat.discr)-1)]
groups <- dat.discr[,ncol(dat.discr)]
results <- fs.stability(vars, 
                        groups, 
                        method = c("plsda", "rf"), 
                        f = 20, 
                        k = 3, 
                        k.folds = 10, 
                        verbose = FALSE)


###################################################
### code chunk number 5: OmicsMarkeR.Rnw:146-147
###################################################
performance.metrics(results)


