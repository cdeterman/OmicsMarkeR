dat.discr <- create.discr.matrix(
    create.corr.matrix(
        create.random.matrix(nvar = 50, 
                             nsamp = 100, 
                             st.dev = 1, 
                             perturb = 0.2)),
    D = 10
)

vars <- dat.discr$discr.mat
groups <- dat.discr$classes

fits <- fs.stability(vars, 
                     groups, 
                     method = c("plsda", "rf"), 
                     f = 10, 
                     k = 3, 
                     k.folds = 10, 
                     verbose = 'none')

feature.table(fits, "plsda")
