\name{modelTuner}
\alias{modelTuner}
\title{Model Tuner}
\usage{
  modelTuner(trainData, guide, method, inTrain, outTrain,
    lev, savePredictions = FALSE, allowParallel = FALSE,
    verbose = FALSE, theDots = NULL)
}
\arguments{
  \item{trainData}{Data used to fit the model}

  \item{guide}{Output from \code{\link{tune.instructions}}.
  Facilitates the optimization by avoiding redundant model
  fitting.}

  \item{inTrain}{Indicies for cross-validated training
  models}

  \item{outTrain}{Indicies for cross-validated testing
  models}

  \item{lev}{Group levels}

  \item{savePredictions}{Logical argument dictating if
  should save the prediction data.  Default
  \code{savePredictions = FALSE}}

  \item{allowParallel}{Logical argument dictating if
  parallel processing is allowed via foreach package}

  \item{verbose}{Logical argument dictating if should print
  progress}

  \item{theDots}{List of additional arguments provided in
  the initial classification and features selection
  function}
}
\value{
  Returns list of fitted models
}
\description{
  Optimizes each model based upon the parameters provided
  either by the internal \code{\link{denovo.grid}} function
  or by the user.
}
\author{
  Charles E. Determan Jr.
}

