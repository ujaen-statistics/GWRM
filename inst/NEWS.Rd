\name{NEWS}
\title{News for Package 'GWRM'}
\encoding{latin1}
                        
\section{Changes in version 2.1 (2015-09-30)}{
                              
\itemize{
  \item Corrected bug in \code{gw} in relation of optimization methods.
  \item \code{residuals.gw} S3 method includes now the possibility of parallel computing in order to construct the simulated envelope.
  }
}

\section{Changes in version 2.0 (2015-04-07)}{
                              
\itemize{
  \item \code{GWRM.fit} is now \code{gw}. The default algorithm to find ML-estimates is L-BFGS-B after nlm and Nelder-Mead previous optimizations; it provides SE estimates of all the parameters of the model.
  \item \code{gw} with \code{k=1} argument corresponds to the classical Waring regression model.
  \item \code{GWRM.display} is now the S3 method \code{summary.gw}.
  \item \code{GWRM.add} is now included in the new S3 method \code{add1.gw} and may be used through \code{step}.
  \item \code{GWRM.stats} is now the S3 method \code{partvar.gw}.
  \item A \code{residuals.gw} S3 method has been included to validate the fits with the analysis of residuals.
  \item NEWS and CITATION files included.
  }
}