\name{tune.instructions}
\alias{tune.instructions}
\title{Model Optimization Instructions}
\usage{
  tune.instructions(method, grid)
}
\arguments{
  \item{method}{Vector of strings indicating which models
  will be fit}

  \item{grid}{A list of parameters grids to be applied to
  the models}
}
\value{
  \item{modelInfo}{List of the following components}.

  \itemize{ \item{scheme: String dictating which looping
  scheme to apply} \item{loop: Dataframe of parameters to
  loop through for each model} \item{model: Information
  regarding parameters of specific model} \item{constant:
  Names of the 'loop' dataframe components} \item{vary:
  Indication of parameters that vary and can access
  recursively}}
}
\description{
  Provides directions for which parameters to loop over
  during tuning.  This becomes important when certain
  models can access 'lower' parameters without running them
  independently.
}
\author{
  Charles E. Determan Jr.
}

