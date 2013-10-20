\name{create.discr.matrix}
\alias{create.discr.matrix}
\title{Discriminatory Multivariate Data Generator}
\usage{
  create.discr.matrix(V, D = 20, l = 1.5, num.groups = 2)
}
\arguments{
  \item{V}{Numeric matrix}

  \item{D}{Number of discriminatory variables induced.
  Default \code{D = 20}}

  \item{l}{Level of discrimination, higher = greater
  separation.  Default \code{l = 1.5}}

  \item{num.groups}{Number of groups in the dataset}
}
\value{
  Matrix of dimension \code{dim(V)+1} with discriminatory
  variables induced and the .classes added to the end of
  the matrix.
}
\description{
  Generates a matrix of dimensions \code{dim(U)} with
  induced correlations.  \code{D} variables are randomly
  selected as discriminatory.  If \code{num.groups = 2}
  then discrimination is induced by adding and subtracting
  values derived from the level of of discrimination,
  \code{l}, for the classes respectively.  Multi-class
  datasets have a few further levels of randomization.  For
  each variable, a random number of the groups are selected
  as discriminating while the remaining groups are not
  altered.  For each discriminatory group, a unique change
  is provided by randomly assigning addition or subtraction
  of the discrimination factor.  For example, if 3 groups
  are selected and two groups are assigned as addition and
  the third subtraction, the second addition is multiplied
  by its number of replicates.  E.g. (1,1,-1) -> (1,2,-1).
  These values are randomized and then multiplied by the
  respective discrimination factor.  The resulting values
  are then added/subtracted from the respective groups.  A
  noise matrix is applied to the final matrix to perturb
  'perfect' discrimination.
}
\author{
  Charles E. Determan Jr.
}

