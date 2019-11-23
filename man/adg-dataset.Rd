\name{adg}
\alias{adg}
\docType{data}

\title{Dataset related to average daily growth performance and abattoir findings in pigs commercial production.}

      \usage{adg}

\description{The case study dataset is about growth performance and abattoir findings in pigs commercial production in a selected set of 15 Canadian farms collected in March 1987.}


\format{An adapted data frame of the original dataset which  consists of 341 observations of 8 variables and a grouping variable (farm).
  \describe{
    \item{farm}{farm ID.}
    \item{female}{sex of the pig (1=female, 0=castrated). }
    \item{age}{days elapsed from birth to slaughter (days).}
    \item{adg}{average daily weight gain (grams).}
    \item{pneumS}{presence of moderate to severe pneumonia.}
    \item{eggs}{presence of fecal/gastrointestinal nematode eggs at time of slaughter.}
    \item{wormCount}{count of nematodes in small intestine at time of slaughter.}
    \item{livdam}{presence of liver damage (parasite-induced white spots).}
    \item{AR}{presence of atrophic rhinitis.}
    }
 }
 
 \references{
 Kratzer, G., Lewis, F.I., Comin, A., Pittavino, M. and Furrer, R. (2019). "Additive Bayesian Network Modelling with the R Package abn". arXiv preprint arXiv:1911.09006.
 
Dohoo, Ian Robert, Wayne Martin, and Henrik Stryhn. Veterinary epidemiologic research. No. V413 DOHv. Charlottetown, Canada: AVC Incorporated, 2003.}

\keyword{datasets}

