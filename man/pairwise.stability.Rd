\name{pairwise.stability}
\alias{pairwise.stability}
\title{Pairwise Stability Metrics}
\usage{
  pairwise.stability(features, stability.metric, nc)
}
\arguments{
  \item{features}{A matrix of selected features}

  \item{stability.metric}{string indicating the type of
  stability metric. Available options are \code{"jaccard"}
  (Jaccard Index/Tanimoto Distance), \code{"sorensen"}
  (Dice-Sorensen's Index), \code{"ochiai"} (Ochiai's
  Index), \code{"pof"} (Percent of Overlapping Features),
  \code{"kuncheva"} (Kuncheva's Stability Measures),
  \code{"spearman"} (Spearman Rank Correlation), and
  \code{"canberra"} (Canberra Distance)}

  \item{nc}{Number of variables in original dataset}
}
\value{
  A list is returned containing: \item{comparisons}{Matrix
  of pairwise comparisons} \item{overall}{The average of
  all pairwise comparisons}
}
\description{
  Conducts all pairwise comparisons of features selected
  following bootstrapping. Also known as the data
  perturbation ensemble approach
}
\author{
  Charles Determan Jr
}
\references{
  He. Z. & Weichuan Y. (2010) \emph{Stable feature
  selection for biomarker discovery}. Computational Biology
  and Chemistry 34 215-225.
}

