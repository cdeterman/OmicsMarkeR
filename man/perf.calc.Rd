\name{perf.calc}
\alias{perf.calc}
\title{Performance Statistics}
\usage{
  perf.calc(data, lev = NULL, model = NULL)
}
\arguments{
  \item{data}{dataframe of predicted (pred) and observed
  (obs) groups}

  \item{lev}{Group levels}

  \item{model}{String indicating which model was initially
  run}
}
\value{
  Returns confusion matrix and ROC performance statistics
  including Accuracy, Kappa, ROC.AUC, Sensitivity,
  Specificity, Positive Predictive Value, and Negative
  Predictive Value
}
\description{
  Calculates confusion matrix and ROC statistics comparing
  the results of the fitted models to the observed groups.
}
\seealso{
  caret function \code{\link{confusionMatrix}}
}

