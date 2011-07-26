%                               -*- Mode: Rd -*- 
% fitbn.Rd --- 
% Author          : Fraser Lewis
% Created On      : 
% Last Modified By: Fraser Lewis
% Last Modified On: 
% Update Count    : 
% Status          : Unknown, Use with caution!
% 

\name{getmarginal}
\alias{getmarginal}
\encoding{latin1}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{Estimation of marginal posterior density in an additive Bayesian network model}

\description{Estimates the marginal posterior density for a single regression coefficient in an additive Bayesian network model}

\usage{
getmarginal(data.df,dag.m,whichnode=NULL, whichvar="constant", 
                  hyper.params=list("mean"=c(0),"var"=c(1000)),
                        verbose=FALSE,post.x=seq(-0.1,0.1,len=100))
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{data.df}{a data frame containing the data used for learning the network, will be coerced to factors if necessary}
  \item{dag.m}{a matrix defining the network structure, a directed acyclic graph, see details for format. Note that colnames and rownames must be set}
  \item{whichnode}{the name of the node - response variable}
  \item{whichvar}{the name of the explanatory variable in whichnode whose posterior density is required. For the intercept term use ``constant''}
  \item{hyper.params}{prior mean and var in Gaussian prior distribution for regression coefficients on logit scale, see details}
  \item{verbose}{logical, if TRUE then extra information is printed to the terminal} 
  \item{post.x}{the set of points in the domain at which to evaluate the posterior density}
}

\details{
  The procedure \code{getmarginal} estimates the posterior density for a regression coefficient within an additive Bayesian network model and is similar to \code{fitabn} but extracts the posterior density rather than the log marginal likelihood. 
  
  In the network structure definition, dag.m, each row represents a node in the network, and the columns in each row define the parents for that particular node, see the example below for the specific format.

Note that currently all variables must be binary.   

}

\value{A matrix with two columns called x and f, x is the value of the variable and f is the corresponding density value.}


\references{Tierney, L. and Kadane, J. B. (1986). Accurate Approximations for Posterior Moments and Marginal Densities.
Journal of the American Statistical Association, Vol. 81, No. 393, pp. 82-86.


  Further information about \bold{abn} can be found at:\cr
  \url{http://www.vetepi.uzh.ch/}.}

\author{
  Fraser Lewis \email{fraseriain.lewis@uzh.ch}
}

\seealso{\code{\link{fitabn}}}

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
## parameter for node X6 baseline/intercept coefficient
\dontrun{getmarginal(sim10varadd,mymodel.dag,whichnode="X6",
                whichvar="constant",hyper.params=list(
                mean=c(0,0,0,0,0,0,0,0,0,0,0),
                var=c(100,100,100,100,100,100,100,100,100,100,100)),
                verbose=FALSE);
}
## parameter for node X6 coefficient for effect X1
\dontrun{getmarginal(sim10varadd,mymodel.dag,whichnode="X6",
                whichvar="X1",hyper.params=list(
                mean=c(0,0,0,0,0,0,0,0,0,0,0),
                var=c(100,100,100,100,100,100,100,100,100,100,100)),
                verbose=FALSE);
}

}

\keyword{models}
