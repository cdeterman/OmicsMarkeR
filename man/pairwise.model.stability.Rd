\name{pairwise.model.stability}
\alias{pairwise.model.stability}
\title{Pairwise Model Stability Metrics}
\usage{
  pairwise.model.stability(features, stability.metric, m,
    k, nc)
}
\arguments{
  \item{features}{A matrix of selected features}

  \item{stability.metric}{string indicating the type of
  stability metric. Avialable options are \code{"jaccard"}
  (Jaccard Index/Tanimoto Distance), \code{"sorensen"}
  (Dice-Sorensen's Index), \code{"ochiai"} (Ochiai's
  Index), \code{"pof"} (Percent of Overlapping Features),
  \code{"kuncheva"} (Kuncheva's Stability Measures),
  \code{"spearman"} (Spearman Rank Correlation), and
  \code{"canberra"} (Canberra Distance)}

  \item{m}{Number of models fit}

  \item{k}{Number of resamples run}

  \item{nc}{Number of original features}
}
\value{
  A list is returned containing: \item{comparisons}{Matrix
  of pairwise comparisons} \item{overall}{The average of
  all pairwise comparisons}
}
\description{
  Conducts all pairwise comparisons of each model's
  selected features selected following bootstrapping.  Also
  known as the function perturbation ensemble approach
}
\author{
  Charles Determan Jr
}
\references{
  He. Z. & Weichuan Y. (2010) \emph{Stable feature
  selection for biomarker discovery}. Computational Biology
  and Chemistry 34 215-225.
}

