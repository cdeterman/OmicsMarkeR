# Create Multivariate Matrices

# Random Multivariate Matrix

# 50 variables, 100 samples, 1 standard devation, 0.2 noise factor

rand.mat <- create.random.matrix(nvar = 50, 
                                 nsamp = 100, 
                                 st.dev = 1, 
                                 perturb = 0.2)


# Induce correlations in a numeric matrix

# Default settings
# minimum and maximum block sizes (min.block.size = 2, max.block.size = 5)
# default correlation purturbation (k=4)
# see ?create.corr.matrix for citation for methods

corr.mat <- create.corr.matrix(rand.mat)


# Induce Discriminatory Variables

# 10 discriminatory variables (D = 10)
# default discrimination level (l = 1.5)
# default number of groups (num.groups=2)
# default correlation purturbation (k = 4)

dat.discr <- create.discr.matrix(corr.mat, D=10)
