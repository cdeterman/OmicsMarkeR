\name{create.disc.matrix}
\alias{create.disc.matrix}
\title{Discriminatory Multivariate Data Generator}
\usage{
  create.disc.matrix(V, D = 20, l = 1.5, num.groups = 2)
}
\arguments{
  \item{V}{Numeric matrix}

  \item{D}{Number of discriminatory variables induced.
  Default \code{D = 20}}

  \item{l}{Level of discrimination, higher = greater
  separation.  Default \code{l = 1.5}}
}
\value{
  Matrix of dimension \code{dim(V)} with discriminatory
  variables induced.
}
\description{
  Generates a matrix of dimensions \code{dim(U)} with
  induced correlations.  Blocks of variables are randomly
  assigned and correlations are induced.  A noise matrix is
  applied to the final matrix to perturb 'perfect'
  correlations.
}
\author{
  Charles E. Determan Jr.
}

