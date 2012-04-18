%                               -*- Mode: Rd -*- 
% fitbn.Rd --- 
% Author          : Fraser Lewis
% Created On      : 
% Last Modified By: Fraser Lewis
% Last Modified On: 
% Update Count    : 
% Status          : Unknown, Use with caution!
% 

\name{fitbn}
\alias{fitbn}
\encoding{latin1}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{Fit a Bayesian Network to discrete data}

\description{Fits a Bayesian network to discrete data calculating the network score and parameter estimates.}

\usage{
fitbn (data.df, dag.m, prior.obs.per.node=NULL, useK2=FALSE, verbose=FALSE)

}

%- maybe also `usage' for other objects documented here.
\arguments{
  \item{data.df}{a data frame containing the data used for learning the network}
  \item{dag.m}{a matrix defining the network structure, a directed acyclic graph, see details for format}
  \item{prior.obs.per.node}{imaginary database size, see details}
  \item{useK2}{logical, choose either K2 metric or BDeu metric, if FALSE BDeu metric is used, see details}
  \item{verbose}{logical, if TRUE then parameters estimates are printed for each and every parent combination at each node.} 
}

\details{
  The procedure \code{fitbn} fits a Bayesian network to data and estimates the (log) marginal likelihood, commonly referred as the network score. The estimation process comprises of evaluating the analytical expression for the (log) marginal likelihood provided in Heckerman et al (1995). Two different forms of parameter prior are supported, the Bayesian Dirichlet equivalence - uniform - metric, BDeu, and the K2 metric. The former assumes that prior data was observed at each individual node in the network and the weight of this data is equally distributed within the parameter prior at each node. This prior data is typically referred to as the imaginary database size. The BDeu metric is likelihood equivalent. The K2 metric uses a weaker prior of Dir(1,1,...) for each parameter at each node in the network. The K2 metric is not likelihood equivalent. More details can be found in Heckerman et al (1995). If useK2 is TRUE then prior.obs.per.node is ignored. The Bayesian network comprises only binomial or multinomial nodes and all variables will be coerced to factors if required. The verbose option can produce copious output which may need to be re-directed to a file (e.g. using sink()).  
  
  In the network structure definition, dag.m, each row represents a node in the network, and the columns in each row define the parents for that particular node, see the example below for the specific format.  

}


\value{A scalar value, the network score}


\references{Heckerman, D. and Geiger, D. and Chickering, D. M. (1995). Learning Bayesian Networks - The Combination of Knowledge and Statistical-Data, Machine Learning, 20 (3), 197-243.

Lewis, F. I., Brulisauer, F. and Gunn, G. J. (2011). Structure discovery in Bayesian networks: An analytical tool for analysing complex animal health data Preventive Veterinary Medicine,  100, 109-115.

  Further information about \bold{abn} can be found at:\cr
  \url{http://www.vetepi.uzh.ch/research/bgm.html}
}

\author{
  Fraser Lewis \email{fraseriain.lewis@uzh.ch}
}


\examples{

## fit DAG to categorical data using K2 metric 
bin.nodes<-c(1,3,4,6,9,10,11,12,15,18,19,20,21,26,27,28,32); 
var33.cat<-var33[,bin.nodes];#categorical nodes only
mydag<-matrix(c(
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v1
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v3
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v4  
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v6  
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v9  
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v10  
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v11  
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v12  
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v15  
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v18 
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v19
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v20 
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v21 
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v26 
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v27 
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v28 
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0  #v32 
              ),byrow=TRUE,ncol=17); 
colnames(mydag)<-rownames(mydag)<-names(var33.cat);#set names
## now fit the model defined in mydag - full independence model
fitbn (data.df=var33.cat, dag.m=mydag,useK2=TRUE);
# this returns the network score goodness of fit = log marginal likelihood

}

\keyword{models}
