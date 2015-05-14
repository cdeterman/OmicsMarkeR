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

# binary class feature ranking
svmrfeFeatureRanking(x = vars,
                     y = groups, 
                     c = 0.1,
                     perc.rem = 10)
