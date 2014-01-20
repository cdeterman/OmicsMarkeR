\name{RPT}
\alias{RPT}
\title{Robustness-Performance Trade-Off}
\usage{
RPT(stability, performance, beta = 1)
}
\arguments{
  \item{stability}{Stability metric i.e. result from
  jaccard, sorensen, etc.}

  \item{performance}{Model performance e.g. accuracy}

  \item{beta}{Relative of importance of stability versus
  performance. Default \code{beta = 1} treats stability and
  performance equally.}
}
\value{
Harmonic mean of robustness and classification performance
}
\description{
A variation on the F-measure (precision and recall) to
assess robustness versus classification performance.
}
\references{
Saeys Y., Abeel T., et. al. (2008) \emph{Machine Learning
and Knowledge Discovery in Databases}. 313-325.
http://link.springer.com/chapter/10.1007/978-3-540-87481-2_21
}

