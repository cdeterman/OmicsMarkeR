dat.discr <- create.discr.matrix(
    create.corr.matrix(
        create.random.matrix(nvar = 50, 
                             nsamp = 100, 
                             st.dev = 1, 
                             perturb = 0.2)),
    D = 10
)$discr.mat


vars <- dat.discr[,1:(ncol(dat.discr)-1)]
groups <- as.factor(dat.discr[,ncol(dat.discr)])

fit <- fit.only.model(X=vars, 
                      Y=groups, 
                      method="plsda", 
                      p = 0.9)
