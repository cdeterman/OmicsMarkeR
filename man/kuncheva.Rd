\name{kuncheva}
\alias{kuncheva}
\title{Kuncheva's Index}
\usage{
  kuncheva(x, y, num.features)
}
\arguments{
  \item{x}{vector of feature names}

  \item{y}{vector of feature names}
}
\value{
  Returns the Kuncheva Index for the two vectors. It takes
  values in [0,1], with 0 meaning no overlap between two
  sets and 1 meaning two sets are identical.
}
\description{
  Calculates Kuncheva's index between two vectors of
  features.  In brief, the closer to 1 the more similar the
  vectors.  The two vectors must have the same cardinality
  (i.e. same length).
}
\references{
  Kuncheva L. (2007) \emph{A stability index for feature
  selection}. Proceedings of the 25th IASTED International
  Multi-Conference: Artificial Intelligence and
  Applications. pp. 390-395.

  He. Z. & Weichuan Y. (2010) \emph{Stable feature
  selection for biomarker discovery}. Computational Biology
  and Chemistry 34 215-225.
}
\seealso{
  \code{\link{kuncheva}}, \code{\link{sorensen}},
  \code{\link{ochiai}}, \code{\link{pof}},
  \code{\link{pairwise.stability}},
  \code{\link{pairwise.model.stability}}
}

