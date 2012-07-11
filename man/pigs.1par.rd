

\name{pigs.1par}
\alias{pigs.1par}
\non_function{}
\title{simulated dataset from a DAG comprising of 13 variables}
\description{
 9011 observations simulated from a DAG with 10 binary variables and 3 continuous variables. A DAG of this data features in the vignette.  
}
\format{A data frame with a mixture of discrete variables each of which is set as a factor and continuous variables.
  \describe{
    \item{D1}{Binary,  conditionally dependent upon D2 }
    \item{D2}{Binary,  conditionally dependent upon D3. }
    \item{D3}{Binary,  conditionally dependent upon D4. }
    \item{D4}{Binary,  conditionally dependent upon Loc.x. }
    \item{D5}{Binary,  conditionally dependent upon D6. }
    \item{D6}{Binary,  conditionally dependent upon D4. }
    \item{D7}{Binary,  conditionally dependent upon Year. }
    \item{D8}{Binary,  conditionally dependent upon D10. }
    \item{D9}{Binary,  conditionally dependent upon D2. }
    \item{D10}{Binary, conditionally dependent upon D9. }
    \item{Year}{Gaussian,  independent. }
    \item{Loc.x}{Gaussian, conditionally dependent upon D7.}
    \item{Loc.y}{Gaussian, conditionally dependent upon D7.}
   
    }
 }

\keyword{datasets}

