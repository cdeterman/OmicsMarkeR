\name{canberra}
\alias{canberra}
\title{Canberra Distance}
\usage{
  canberra(x, y)
}
\arguments{
  \item{x}{numeric vector of ranks}

  \item{y}{numeric vector of ranks with compatible length
  to x}
}
\value{
  Returns the canberra distance for the two vectors
}
\description{
  Calculates canberra distance between two vectors.  In
  brief, the higher the canberra distance the greater the
  'distance' between the two vectors (i.e. they are less
  similar).
}
\references{
  Jurman G., Merler S., Barla A., Paoli S., Galea A., &
  Furlanello C. (2008) \emph{Algebraic stability indicators
  for ranked lists in molecular profiling}. Bioinformatics
  24(2): 258-264.

  He. Z. & Weichuan Y. (2010) \emph{Stable feature
  selection for biomarker discovery}. Computational Biology
  and Chemistry 34 215-225.
}

