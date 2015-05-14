\name{News}
\title{News for Package 'OmicsMarkeR'}

\section{Version 1.1.0}{
    \itemize{
        \item Changes:
            \itemize{
                \item the \code{modelList} has been added to provide the user
                a list of currently implemented methods.
                \item The 'verbose' argument in fs.stability, 
                fs.ensembl.stability, & fit.only.model has been changed
                to a character option indicating the extent of verbose output.
                \item Canberra stability has been added as the previous
                implementation was not compatible with RPT. See 
                \code{?canberra.stability} for more details.
                \item Many more tests have been implemented to make the
                package more stable.
            }
    }
}

\section{Initial release Version 0.99.3}{
    \itemize{
        \item FEATURES IN PROGRESS:
        \itemize{
            \item Access to fitted models (averaged or all bootstrapped 
            results?)
            \item Easy graphics access (scores/loadings plots, variable 
            importance plots, etc.)
            \item Summary graphics (across models)
            \item Database searching (HMDB, MMCD, Metlin, LipidMaps, etc.)
            \item Additional algorithms (neural nets, bayesian, etc.)
            \item Additional ensemble methods (bayesian, boosting, etc.)
        }
    }
}
    

    
