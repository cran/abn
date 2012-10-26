% abn-internal.Rd --- 
% Author           : Fraser Lewis
% Last modified on : 03/08/2012 

\name{calc.node.mlik.inla}

\alias{calc.node.mlik.inla}

\title{internal function for individual node computation via INLA}

\description{
  This functions is for internal use only and not meant to be called by the user.
}

\usage{
calc.node.mlik.inla(child.loc,dag.m.loc,data.df.loc,data.dists.loc,ntrials.loc,exposure.loc,compute.fixed.loc,
                         mean.intercept.loc,prec.intercept.loc,mean.loc,prec.loc,loggam.shape.loc,loggam.inv.scale.loc,verbose.loc)
}
\arguments{ 
  \item{child.loc}{row index of current node in dag.m}
  \item{dag.m.loc}{a matrix defining the network structure, a directed acyclic graph (DAG), see details for format. Note that colnames and rownames must be set}
  \item{data.df.loc}{a data frame containing the data used for learning the network, binary variables must be declared as factors}
  \item{data.dists.loc}{a named list giving the distribution for each node in the network, see details}
  \item{ntrials.loc}{a numeric vector giving the total number of trials, only applicable if the data comprise one or more binary variables. see details}
  \item{exposure.loc}{a numeric vector, giving the unit of exposure, only applicable if the data comprise one or more poisson variables, see details}
  \item{compute.fixed.loc}{a logical flag, set to TRUE if also want the marginal posterior distributions for the model parameters as well as the marginal likelihood}
  \item{mean.intercept.loc}{the prior mean of the Gaussian intercept term for each node}
  \item{prec.intercept.loc}{the prior precision of the Gaussian intercept term for each node}
  \item{mean.loc}{the prior mean in the Gaussian prior for each covariate at each node}
  \item{prec.loc}{the prior precision in the Gaussian prior for each covariate at each node}
  \item{loggam.shape.loc}{the prior shape parameter in the Gamma prior for the precision parameters at each node}
  \item{loggam.inv.scale.loc}{the prior inverse scale (=1/scale) parameter in the Gamma prior for the precision parameters at each node}
  \item{verbose.loc}{if true then provides some additional output, in particular the call used to inla()}
}


\value{A named list. If compute.fixed=FALSE then each entry in the list is the log marginal likelihood for each node, else each entry in this list is itself a list comprising of all the different outputs provided by \code{inla}.}

\details{
This function is called in \code{\link{fitabn}} and \code{\link{buildscorecache}} and deals with the actual call to INLA which builds up an appropriate string and is then passed into eval(parse()). The arguments all have the .loc postix due to a scoping issue using eval(parse()) which does not seem to use the parent.frame/calling environment scope. 

}

\author{
  Fraser Lewis \email{fraseriain.lewis@uzh.ch}.
}


\keyword{internal}
