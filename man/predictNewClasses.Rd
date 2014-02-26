\name{predictNewClasses}
\alias{predictNewClasses}
\title{Class Prediction}
\usage{
predictNewClasses(modelFit, method, orig.data, newdata, parms = NULL)
}
\arguments{
  \item{method}{String of the model to be evaluated}

  \item{modelFit}{The fitted model being evaluated}

  \item{orig.data}{The orginal data before subsetting
  training sets.  Required to have the 'observed' group
  membership}

  \item{indicies}{The indicies for the training subsets}

  \item{newdata}{The testing data to predict group
  membership}

  \item{parms}{Optional alternate parameters being fit to
  the model}
}
\value{
Returns a list of predicted group membership
}
\description{
This function evaluates a single fitted model and returns
the predicted group memberships of new data.
}

