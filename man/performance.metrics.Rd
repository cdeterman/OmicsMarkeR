\name{performance.metrics}
\alias{performance.metrics}
\title{Performance Metrics of fs.stability or fs.ensembl.stability object}
\usage{
performance.metrics(fit.model, digits = max(3, getOption("digits") - 3))
}
\arguments{
  \item{fit.model}{An fs.stability or fs.ensembl.stability
  object}

  \item{digits}{How many digits to round values}
}
\value{
Dataframe of performance statistics by model
}
\description{
This will provide a concise data.frame of confusion matrix
and ROC statistics from the results of \code{fs.stability}
or \code{fs.ensembl.stability}.
}
\author{
Charles E. Determan Jr.
}

