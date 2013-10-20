\name{tune}
\alias{tune}
\title{Model Tuning and Metrics}
\usage{
  tune(trainVars, trainGroup, method, k.folds = 10,
    repeats = 3, res = 3, grid = NULL, metric = "Accuracy",
    savePerformanceMetrics = NULL, allowParallel = FALSE,
    verbose = FALSE, theDots = NULL)
}
\arguments{
  \item{trainVars}{Data used to fit the model}

  \item{trainGroup}{Group identifiers for the training
  data}

  \item{method}{A vector of strings listing models to be
  optimized}

  \item{k.folds}{Number of folds generated during
  cross-validation.  Default \code{"k.folds = 10"}}

  \item{repeats}{Number of times cross-validation repeated.
  Default \code{"repeats = 3"}}

  \item{res}{Resolution of model optimization grid.
  Default \code{"res = 3"}}

  \item{grid}{Optional list of grids containing parameters
  to optimize for each algorithm. Default \code{"grid =
  NULL"} lets function create grid determined by
  \code{"res"}}

  \item{metric}{Criteria for model optimization.  Available
  options are \code{"Accuracy"} (Predication Accuracy),
  \code{"Kappa"} (Kappa Statistic), and \code{"AUC-ROC"}
  (Area Under the Curve - Receiver Operator Curve)}

  \item{savePerformanceMetrics}{Logical argument dictating
  if should save the prediction data.  Default
  \code{savePerformanceMetrics = FALSE}}

  \item{allowParallel}{Logical argument dictating if
  parallel processing is allowed via foreach package}

  \item{verbose}{Logical argument if should output
  progress}

  \item{theDots}{List of additional arguments provided in
  the initial classification and features selection
  function}
}
\value{
  Basically a list with the following elements:

  \item{method}{Vector of strings listing models that were
  optimized}

  \item{performance}{Performance generated internally to
  optimize model}

  \item{bestTune}{List of paramaters chosen for each model}

  \item{dots}{List of extra arguments initially provided}

  \item{metric}{Criteria that was used for model
  optimization}

  \item{finalModels}{The fitted models with the 'optimum'
  parameters}

  \item{performance.metrics}{The performance metrics
  calculated internally for each resulting prediction}

  \item{tune.metrics}{The results from each tune}

  \item{perfNames}{The names of the performance metrics}

  \item{comp.catch}{If the optimal PLSDA model contains
  only 1 component, the model must be refit with 2
  components.  This catches the 1 component parameter so
  feature selection and further performance analysis can be
  conducted on the 1 component.}
}
\description{
  Optimizes each model based upon the parameters provided
  either by the internal \code{\link{denovo.grid}} function
  or by the user.
}
\author{
  Charles E. Determan Jr.
}

