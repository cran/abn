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

\title{Fit an additive Bayesian network model to binary data}

\description{Fits an additive Bayesian network to binary data using a logistic link function and is equivalent to a (Bayesian) multi-dimensional logistic regression model.}

\usage{
fitabn (data.df,dag.m,hyper.params=list("mean"=c(0),"var"=c(1000)),
               verbose=FALSE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{data.df}{a data frame containing the data used for learning the network, will be coerced to factors if necessary}
  \item{dag.m}{a matrix defining the network structure, a directed acyclic graph, see details for format. Note that colnames and rownames must be set}
  \item{hyper.params}{prior mean and var in Gaussian prior distribution for regression coefficients on logit scale, see details}
  \item{verbose}{logical, if TRUE then extra information is printed to the terminal} 
}

\details{
  The procedure \code{fitabn} fits an additive Bayesian network model to data using a logistic link function, estimating the log marginal likelihood for the model using the laplace method for approximating intergals (see references). This model is equivalent to a (Bayesian) multi-dimensional logistic regression model. This is designed to complement the usual (conjugate) contingency table Bayesian network models (which are provided in \code{fitbn}). An additive model uses far fewer parameters but is computationally much more demanding due to the lack of conjugacy, and hence an analytical solution to the marginal likelihood. Currently the parameter priors are assumed independent and do not vary with the node (response variable) in the network, and it only makes sense for these to be uninformative e.g. of zero mean and high variance. To override the default settings of means of zero and variance of 1000 then hyper.params must be passed a list with two named vectors, mean and var, where these must contain d+1 entries, where d is the number of columns in dag.m, and the first term in each vector is the prior for the constant term which is implicitly present in all nodes and is the baseline (corner point constraints). The marginal posterior distributions for each parameter are available through \code{getmarginal}. 
  
  In the network structure definition, dag.m, each row represents a node in the network, and the columns in each row define the parents for that particular node, see the example below for the specific format.

Note that currently all variables must be binary.   

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
data(sim10varadd);
##  a matrix of arcs to be banned from the search process
# - do not ban any arcs - empty matrix
mymodel.dag<-matrix(c( 0,0,0,0,0,0,0,0,0,0, ## X1
                       0,0,0,0,1,0,0,0,0,0, ## X2|X5
                       0,0,0,0,0,0,0,0,0,0, ## X3
                       0,0,0,0,0,0,0,0,0,0, ## X4
                       0,0,0,0,0,0,0,0,0,0, ## X5
                       1,0,0,0,0,0,1,0,1,0, ## X6|X1:X7:X9
                       0,0,1,0,0,0,0,0,0,0, ## X7|X3
                       0,0,0,0,0,0,0,0,0,0, ## X8
                       0,0,0,0,0,0,0,0,0,0, ## X9
                       0,0,0,0,0,0,1,0,0,0),## X10|X7
                       byrow=TRUE,ncol=10);

colnames(mymodel.dag)<-names(sim10varadd);# must set names
rownames(mymodel.dag)<-names(sim10varadd);# must set names

\dontrun{fitabn(sim10varadd,mymodel.dag,hyper.params=list(
                mean=c(0,0,0,0,0,0,0,0,0,0,0),
                var=c(100,100,100,100,100,100,100,100,100,100,100)),
                verbose=FALSE);
}



}

\keyword{models}
