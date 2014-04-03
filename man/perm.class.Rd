\name{perm.class}
\alias{perm.class}
\title{Monte Carlo Permutation of Model Performance}
\usage{
perm.class(fs.model = NULL, X, Y, method, k.folds = 5,
  metric = "Accuracy", nperm = 10, allowParallel = FALSE, ...)
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

  \item{k.folds}{How many and what fractions of dataset
  held-out for prediction (i.e. 3 = 1/3, 10 = 1/10, etc.)}

  \item{metric}{Performance metric to assess.  Available
  options are \code{"Accuracy"}, \code{"Kappa"}, and
  \code{"ROC.AUC"}.}

  \item{nperm}{Number of permutations, default \code{nperm
  = 10}}

  \item{allowParallel}{Logical argument dictating if
  parallel processing is allowed via foreach package.
  Default \code{allowParallel = FALSE}}

  \item{...}{Extra arguments that the user would like to
  apply to the models}
}
\value{
\item{p.value}{Resulting p-value of permuation test}
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
Guo Y., et. al. (2010) \emph{Sample size and statistical
power considerations in high-dimensionality data settings:
a comparative study of classification algorithms}. BMC
Bioinformatics 11:447.
}

