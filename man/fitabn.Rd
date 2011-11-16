%                               -*- Mode: Rd -*- 
% fitbn.Rd --- 
% Author          : Fraser Lewis
% Created On      : 
% Last Modified By: Fraser Lewis
% Last Modified On: 
% Update Count    : 
% Status          : Unknown, Use with caution!
% 

\name{fitabn}
\alias{fitabn}
\encoding{latin1}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{Fit an additive Bayesian network model to continuous or binary data}

\description{Fits an additive Bayesian network to continuous or binary data and is equivalent to a (Bayesian) multi-dimensional regression model. A logistic link is used for binary variables.}

\usage{
fitabn (data.df,dag.m,hyper.params=list("mean"=c(0),
                                        "sd"=c(sqrt(1000)),
                                        "shape"=c(0.001),
                                        "scale"=c(1/0.001)),
        max.iters=100,epsabs=1e-7, verbose=FALSE,error.verbose=FALSE)
}
    
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{data.df}{a data frame containing the data used for learning the network, will be coerced to factors if necessary}
  \item{dag.m}{a matrix defining the network structure, a directed acyclic graph, see details for format. Note that colnames and rownames must be set}
  \item{hyper.params}{prior means and standard deviations in Gaussian prior distribution for regression coefficients, see details}
  \item{max.iters}{maximum number of iterations used to find root in Laplace approximation}
  \item{epsabs}{minimum absolute error of root estimate in Laplace approximation}
  \item{verbose}{logical, if TRUE then extra information is printed to the terminal}
 \item{error.verbose}{logical, if TRUE then extra information is printed to the terminal about any numerical convergence issues in the Laplace approximation}
}

\details{
 The procedure \code{fitabn} fits an additive Bayesian network model to data where each node (variable) can be either binary (binomial) or continuous (Gaussian) and the model comprises of a set of conditionally independent generalized linear regressions. For a variable to be treated as binary it has to be a factor (with two levels) otherwise it will be assumed to be Gaussian distributed. For binary nodes a logistic link function is used. Estimating the log marginal likelihood for the model uses the Laplace method for approximating intergals (see references). This additive network model is designed to complement the usual (conjugate) contingency table Bayesian network models (which are provided in \code{fitbn}). An additive model uses far fewer parameters but is computationally much more demanding due to the lack of conjugacy, and hence an analytical solution to the marginal likelihood. Independent Gaussian priors are assumed for all coefficients and independent Gamma priors for the precision (=1/variance) terms in each Gaussian node. The shape and scale parameters as are defined in \code{rgamma}. To override the default hyper-parameters then hyper.params must be passed a list with sufficient entries in each vector element. For example the named vectors, mean and sd, must contain d+1 entries, where d is the number of columns in dag.m, where the first term in each vector is the prior for the (node specific) constant/intercept term which is implicitly present in all nodes and is the baseline (using corner point constraints). The length of shape and scale vectors must match the number of Gaussian nodes i.e. the number of variables in data.df which are not factors. If the data comprise of all binary nodes then shape and scale should be passed as NULL. Marginal posterior distributions for each parameter are available through \code{getmarginal}. 
  
In the network structure definition, dag.m, each row represents a node in the network, and the columns in each row define the parents for that particular node, see the example below for the specific format.

}

\value{A scalar value, the network score - log marginal likelihood}


\references{Tierney, L. and Kadane, J. B. (1986). Accurate Approximations for Posterior Moments and Marginal Densities.
Journal of the American Statistical Association, Vol. 81, No. 393, pp. 82-86.

Rijmen, F. (2008). Bayesian networks with a logistic regression model for the conditional probabilities. International Journal of Approximate Reasoning, Vol 48, No. 2, pp. 656-666

  Further information about \bold{abn} can be found at:\cr
  \url{http://www.vetepi.uzh.ch/}}

\author{
  Fraser Lewis \email{fraseriain.lewis@uzh.ch}
}

\seealso{\code{\link{fitbn}}}

\examples{
cts.nodes<-c(2,5,7,8,13,14,16,17,22,23,24,25,29,30,31,33);
var33.cts<-var33[,cts.nodes];#drop categorical nodes
mydag<-matrix(c(
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v2
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v5  
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v7  
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v8 
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v13  
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v14  
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v16  
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v17  
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v22 
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v23
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v24 
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v25 
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v29 
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v30 
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v31 
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0  #v33 
              ),byrow=TRUE,ncol=16); 
colnames(mydag)<-rownames(mydag)<-names(var33.cts);#set names
## now fit the model defined in mydag - full independence
fitabn (data.df=var33.cts,dag.m=mydag,verbose=TRUE);
## uses default priors of N(mu=0,var=1000), 1/var=Gamma(0.001,1/0.001)
# this is the network score goodness of fit = log marginal likelihood


}

\keyword{models}
