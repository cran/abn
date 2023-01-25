%                               -*- Mode: Rd -*- 
% abn-package.Rd --- 
% Author          : Gilles Kratzer
% Created On      : 22/05/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{. abn .}
\alias{overview}
\alias{abn}
\alias{abn-package}
\docType{package}
\title{abn Package}
\description{\code{abn} is a collection of functions for fitting, selecting/learning, analysing, reporting Additive Bayesian Networks.
}
\section{General overview}{What is \pkg{abn}:\cr

Bayesian network modeling is a data analysis technique that is ideally suited to messy, highly correlated, and complex datasets. This methodology is somewhat distinct from other forms of statistical modeling in that its focus is on structure discovery - determining an optimal graphical model that describes the inter-relationships in the underlying processes which generated the data. It is a multivariate technique and can used for one or many dependent variables. This is a data-driven approach, as opposed to, rely only on subjective expert opinion to determine how variables of interest are inter-related (for example, structural equation modeling).

The R package abn is designed to fit additive Bayesian models to observational datasets. It contains routines to score Bayesian Networks based on Bayesian or information-theoretic formulation of generalized linear models. It is equipped with exact search and greedy search algorithms to select the best network. The Bayesian implementation supports random effects to control for one layer clustering. It supports a possible mixture of continuous, discrete, and count data and input of prior knowledge at a structural level.

The R package abn requires the R package Rgraphviz to work well. It is  store outside of CRAN; see \sQuote{Examples} for the code to install the last version.

The web pages \url{http://r-bayesian-networks.org} provide futher case
studies. See also the files provided in the package directories
\code{inst/bootstrapping_example} and \code{inst/old_vignette} for more details.    

}
\references{
Kratzer, G., Lewis, F.I., Comin, A., Pittavino, M. and Furrer, R. (2023). Additive Bayesian Network Modelling with the R Package abn. 
\emph{Journal of Statistical Software}, 105(8), 1--41, \doi{10.18637/jss.v105.i08}.

Lewis, F. I., and Ward, M. P. (2013). "Improving epidemiologic data analyses through multivariate regression modeling". 
\emph{Emerging Themes in Epidemiology}, 10(1), 4.

}


\examples{
## Citations:
print(citation('abn'), bibtex=TRUE)


## Installing the R package Rgraphviz:
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("Rgraphviz")

## README.md in the directory `bootstrapping_example/`:
# edit(file=paste0( path.package('abn'),'/bootstrapping_example/README.md'))
}
\author{Fraser Iain Lewis and Gilles Kratzer}
\keyword{documentation}
\keyword{package}
