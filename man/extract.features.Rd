\name{extract.features}
\alias{extract.features}
\title{Feature Extraction}
\usage{
  extract.features(x, dat = NULL, grp = NULL, method,
    model.features = FALSE, bestTune = NULL, f,
    comp.catch = NULL)
}
\arguments{
  \item{x}{Previously fitted model}

  \item{dat}{Numeric variable data used for fitted models
  (In appropriate format)}

  \item{method}{String indicating the INDIVIDUAL model
  being extracted from}

  \item{model.features}{Logical argument dictating if
  features selected determined by models instead of user
  determined number of features.}

  \item{bestTune}{If \code{model.features = TRUE}, must
  provide the parameter at which to extract features from
  the model.}

  \item{comp.catch}{An internal check for plsda models.  If
  the optimal model contains only 1 component, the ncomp
  paramter must be set to 2 for the model.  However,
  features are still extracted only from the first
  component.}
}
\value{
  Returns list of the features selected from the fitted
  model.
}
\description{
  Extracts features from models that have been previously
  fit.
}

