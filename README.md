# `"OmicsMarkeR"`

[![Travis-CI Build Status](https://travis-ci.org/cdeterman/OmicsMarkeR.png?branch=master)](https://travis-ci.org/cdeterman/OmicsMarkeR) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/cdeterman/OmicsMarkeR?branch=master)](https://ci.appveyor.com/project/cdeterman/OmicsMarkeR)[![Coverage Status](https://coveralls.io/repos/cdeterman/OmicsMarkeR/badge.svg)](https://coveralls.io/r/cdeterman/OmicsMarkeR)

`OmicsMarkeR` is an R package that provides functions for classification and 
feature selection of 'omics' level datasets.


## Motivation

During my studies as a developing Systems Biologist I discovered there were 
often varied techniques to answer the same initial question, how can I classify
high-dimensional data (i.e. metabolomics, proteomics, transcriptomics)?  
A second question usually posed in Biomarker investigations was which features 
are most important to such classification.

I initially pursued the repositories of CRAN and Bioconductor.  I discovered 
such wonderful packages such as `caret` (which I highly recommend); however, 
I was unable to find a means of systematically running multiple algorithms in 
addition to stability metrics to provide confidence with features identified as
important.  This is critical as there seemed little practical benefit to 
classifying 2+ groups if the features identified varied between each test.

In my readings, I came upon an excellent chapter in the Lecture Notes of 
Computer Science Vol. 5212 entitled 'Robust Feature Selection Using Ensemble 
Feature Selection Techniques' by Yvan Sayes, Thomas Abeel, and Yves Van de Peer.
From this chapter I decided to build this package, a tool to provide multiple 
multivariate classification and feature selection techniques complete with 
multiple stability metrics and aggregation techniques.  In this manner, this 
package provides a way to systematically compare both data perturbation and 
function perturbation ensemble techniques complete with a harmonic mean of 
feature robustness and classification performance to evaluate the optimal model 
for the individual dataset.  This following David Wolpert's 'No Free Lunch
Theorem' as there is no single model that is appropriate for all problems.

I have made every effort to cite articles in which either the original 
technique was developed or applied. The interested reader, as well you should 
be, is highly encouraged to seek out these articles.


## Installation

Stable version [Bioconductor](http://www.bioconductor.org/packages/release/bioc/html/OmicsMarkeR.html)
```r
source("http://bioconductor.org/biocLite.R")
biocLite("OmicsMarkeR")
```

## Features in Progress
1. Access to fitted models (averaged or all bootstrapped results?)
2. Easy graphics access (scores/loadings plots, variable importance plots, etc.)
3. Summary graphics (across models)
4. Database searching (HMDB, MMCD, Metlin, LipidMaps, etc.)
5. Additional algorithms
6. Additional ensemble methods (bayesian, boosting, etc.)

