\name{perm.features}
\alias{perm.features}
\title{Feature Selection via Monte Carlo Permutation}
\usage{
perm.features(fs.model = NULL, X, Y, method, sig.level = 0.05, nperm = 10,
  allowParallel = FALSE, ...)
}
\arguments{
  \item{fs.model}{Object containing results from
  \code{fs.stability}}

  \item{X}{A scaled matrix or dataframe containing numeric
  values of each feature}

  \item{Y}{A factor vector containing group membership of
  samples}

  \item{method}{A vector listing models to be fit.
  Available options are \code{"plsda"} (Partial Least
  Squares Discriminant Analysis), \code{"rf"} (Random
  Forest), \code{"gbm"} (Gradient Boosting Machine),
  \code{"svm"} (Support Vector Machines), \code{"glmnet"}
  (Elastic-net Generalized Linear Model), and \code{"pam"}
  (Prediction Analysis of Microarrays)}

  \item{sig.level}{Desired significance level for features,
  default \code{sig.level = .05}}

  \item{nperm}{Number of permutations, default \code{nperm
  = 10}}

  \item{allowParallel}{Logical argument dictating if
  parallel processing is allowed via foreach package.
  Default \code{allowParallel = FALSE}}
}
\value{
\item{sig.level}{User-specified significance level}

\item{num.sig.features}{Number of significant features}

\item{sig.features}{Dataframe of significant features}
}
\description{
Applies Monte Carlo permutations to user specified models.
The user can either use the results from
\code{fs.stability} or provide specified model parameters.
}
\author{
Charles Determan Jr.
}
\references{
Wongravee K., et. al. (2009) \emph{Monte-Carlo methods for
determining optimal number of significant variables.
Application to mouse urinary profiles}. Metabolomics
5:387-406.
}

