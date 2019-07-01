%                               -*- Mode: Rd -*- 
% expit.Rd --- 
% Author          : Gilles Kratzer
% Created On      : 17.07.2018
% Last Modified By: GK (rd file + naming + examples)
% Last Modified On: 
% Update Count    : 
% Status          : Unknown, Use with caution!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{expit}
\alias{expit}
\alias{logit}
\alias{odds}

\encoding{latin1}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{Expit, Logit and odds}

\description{Compute the expit and logit of a numerical vector. Transform odds to probability.}

\usage{
expit(x)
logit(x)
odds(x)
}

\arguments{
  \item{x}{vector of real values.}
}

\details{
\code{logit} computes the logit function:
\deqn{logit(p) = \log\frac{p}{1-p}}{logit(p) = log(p/(1-p)),}
\code{expit} computes the expit function:
\deqn{expit(x) = \frac{e^x}{1+e^x}}{expit(x) = e^x/(1+e^x).}
\code{odds} transform an odd into a probability.
\deqn{odds(x) = \frac{x}{1-x}}{odds(x) = x/(1-x)}
Those functions become numerically unstable if evaluated at the edge or the definition range.
}
\value{
A real vector corresponding to the expit, the logit or the odds of the input values.
}

\author{Gilles Kratzer}

\keyword{manip}
