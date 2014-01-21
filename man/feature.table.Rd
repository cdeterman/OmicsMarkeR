\name{feature.table}
\alias{feature.table}
\title{Feature Consistency Table}
\usage{
feature.table(features, method = NULL)
}
\arguments{
  \item{features}{A fs.stability fitted object}

  \item{method}{Algorithm of interest Available options are
  \code{"plsda"} (Partial Least Squares Discriminant
  Analysis), \code{"rf"} (Random Forest), \code{"gbm"}
  (Gradient Boosting Machine), \code{"svm"} (Support Vector
  Machines), \code{"glmnet"} (Elastic-net Generalized
  Linear Model), and \code{"pam"} (Prediction Analysis of
  Microarrays)}
}
\value{
A data frame containing: \item{features}{Features
identified by model} \item{consistency}{Number of
iterations feature was identified}
\item{frequency}{Frequency of iterations the feature was
identified}
}
\description{
Extracts and sorts the features identified for a given
method.
}
\author{
Charles Determan Jr
}

