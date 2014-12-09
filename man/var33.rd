

\name{var33}
\alias{var33}
\non_function{}
\title{simulated dataset from a DAG comprising of 33 variables}
\description{
 250 observations simulated from a DAG with 17 binary variables and 16 continuous. A DAG of this data features in the vignette. Note that the conditional dependence relations given are those in the population and may differ in the realization of 250 observations. 
}
\format{A data frame with a mixture of discrete variables each of which is set as a factor and continuous variables.
  \describe{
    \item{v1}{Binary, independent. }
    \item{v2}{Gaussian, conditionally dependent upon v1. }
    \item{v3}{Binary, independent. }
    \item{v4}{Binary, conditionally dependent upon v3. }
    \item{v5}{Gaussian, conditionally dependent upon v6. }
    \item{v6}{Binary, conditionally dependent upon v4 and v7. }
    \item{v7}{Gaussian,  conditionally dependent upon v8. }
    \item{v8}{Gaussian,  conditionally dependent upon v9. }
    \item{v9}{Binary, conditionally dependent upon v10. }
    \item{v10}{Binary, independent. }
    \item{v11}{Binary,  conditionally dependent upon v10, v12 and v19. }
    \item{v12}{Binary, independent.}
    \item{v13}{Gaussian, independent.}
    \item{v14}{Gaussian,  conditionally dependent upon v13. }
    \item{v15}{Binary, conditionally dependent upon v14 and v21. }
    \item{v16}{Gaussian, independent. }
    \item{v17}{Gaussian, conditionally dependent upon v16 and v20. }
    \item{v18}{Binary,  conditionally dependent upon v20. }
    \item{v19}{Binary,  conditionally dependent upon v20. }
    \item{v20}{Binary, independent. }
    \item{v21}{Binary, conditionally dependent upon v20. }
    \item{v22}{Gaussian, conditionally dependent upon v21. }
    \item{v23}{Gaussian, conditionally dependent upon v21. }
    \item{v24}{Gaussian, conditionally dependent upon v23. }
    \item{v25}{Gaussian, conditionally dependent upon v23 and v26. }
    \item{v26}{Binary, conditionally dependent upon v20. }
    \item{v27}{Binary, independent. }
    \item{v28}{Binary, conditionally dependent upon v27, v29 and v31. }
    \item{v29}{Gaussian, independent. }
    \item{v30}{Gaussian,  conditionally dependent upon v29. }
    \item{v31}{Gaussian, independent. }
    \item{v32}{Binary,  conditionally dependent upon v21, v29 and v31. }
    \item{v33}{Gaussian,  conditionally dependent upon v31. }
    }
 }

\keyword{datasets}

