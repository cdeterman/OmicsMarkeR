\name{params}
\alias{params}
\title{Model Parameters and Properties}
\usage{
params(method = NULL)
}
\arguments{
  \item{method}{A vector of strings listing the models to
  be returned}
}
\value{
Returns a dataframe of the following components:

method A vector of strings listing models returned

parameter A vector of possible parameters to be optimized

label A vector of the names for each possible parameter

seq A logical indicator if the parameter is sequential in
the model (i.e. if model is able to fit all 'lower'
parameters simultaneously)
}
\description{
Provides a list of the models with their respective
parameters and properties.
}

