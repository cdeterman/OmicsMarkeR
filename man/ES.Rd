\name{ES}
\alias{ES}
\title{Ensemble Stability Aggregation}
\usage{
ES(efs, f)
}
\arguments{
  \item{efs}{A matrix of selected features}

  \item{f}{The number of features desired.  If rank
  correlation desired, f = NULL}
}
\value{
\item{agg}{Aggregated list of features}
}
\description{
Compiles matrix of ranked features via ensemble stability
aggregation
}
\author{
Charles Determan Jr
}
\references{
Meinshausen N., Buhlmann P. (2010) \emph{Stability
selection}. J.R. Statist. Soc. B. 72:4 417-473.
}

