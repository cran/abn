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

Bayesian network modelling is a data analysis technique which is ideally suited to messy, highly correlated and complex datasets. This methodology is rather distinct from other forms of statistical modelling in that its focus is on structure discovery - determining an optimal graphical model which describes the inter-relationships in the underlying processes which generated the data. It is a multivariate technique and can used for one or many dependent variables. This is a data driven approach, as opposed to, rely only on subjective expert opinion to determine how variables of interest are inter-related (for example: structural equation modelling).

The R package abn is designed to fit additive Bayesian models to observational datasets. It contains routines to score Bayesian Networks based on Bayesian or information-theoretic formulation of generalized linear models. It is equipped with exact search and greedy search algorithms to select the best network. The Bayesian implementation supports random effects to control for one layer clustering. It supports a possibly mixture of continuous, discrete and count data and input of prior knowledge at a structural level.
}
\references{
  Lewis, F. I., and Ward, M. P. (2013). Improving epidemiologic data analyses through multivariate regression modelling. Emerging themes in epidemiology,
  10(1), 4.
  
  Kratzer, G., Pittavino, M, Lewis, F. I., and Furrer, R., (2017). abn: an R package for modelling multivariate data using additive Bayesian networks. R
  package version 1.1.  https://CRAN.R-project.org/package=abn
  }


\examples{
## Citations:
citation('abn')
}
\author{Fraser Iain Lewis and Gilles Kratzer}
\keyword{documentation}
\keyword{package}
