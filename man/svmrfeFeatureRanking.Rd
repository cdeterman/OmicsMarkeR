\name{svmrfeFeatureRanking}
\alias{svmrfeFeatureRanking}
\title{SVM Recursive Feature Extraction (Binary)}
\usage{
  svmrfeFeatureRanking(x, y, c, perc.rem = 10)
}
\arguments{
  \item{x}{A matrix where each column represents a feature
  and each row represents a sample}

  \item{y}{A vector of labels corresponding to each
  sample's group membership}

  \item{c}{A numeric value corresponding to the 'cost'
  applied during the svm model fitting. This can be
  selected by the user if using this function directly or
  is done internally.}

  \item{perc.rem}{A numeric value indicating the percent of
  features removed during each iteration. Default
  \code{perc.rem = 10}.}
}
\value{
  Vector of features ranked from most important to least
  important.
}
\description{
  This conducts feature selection for Support Vector
  Machines models via recursive feature extraction.  This
  returns a vector of the features in x ordered by
  relevance. The first item of the vector has the index of
  the feature which is more relevant to perform the
  classification and the last item of the vector has the
  feature which is less relevant. This function is specific
  to Binary classification problems,
}
\references{
  Guyon I. et. al. (2010) \emph{Gene Selection for Cancer
  Classification using Support Vector Machines}. Machine
  Learning 46 389-422.
}
\seealso{
  \code{\link{svmrfeFeatureRankingForMulticlass}}
}

