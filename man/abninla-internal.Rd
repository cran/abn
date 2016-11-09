% abn-internal.Rd --- 
% Author           : Fraser Lewis
% Maintainer       : Marta Pittavino
% Last modified on : 26/09/2014 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{abninla-internal}

\alias{calc.node.inla.glm}
\alias{calc.node.inla.glmm}

\title{abn internal functions}

\description{
These are functions for either internal use or error checking and are
not meant to be called by the user.
}

\usage{

calc.node.inla.glm(child.loc,dag.m.loc,data.df.loc,data.dists.loc,
                   ntrials.loc,exposure.loc,compute.fixed.loc,
                   mean.intercept.loc,prec.intercept.loc,mean.loc,prec.loc,
                   loggam.shape.loc,loggam.inv.scale.loc,verbose.loc)
calc.node.inla.glmm(child.loc,dag.m.loc,data.df.loc,data.dists.loc,
                    ntrials.loc,exposure.loc,compute.fixed.loc,
                    mean.intercept.loc,prec.intercept.loc,mean.loc,prec.loc,
                    loggam.shape.loc,loggam.inv.scale.loc,verbose.loc)
          }

\details{
calc.node.inla.glm() and calc.node.inla.glmm() are internal wrappers to INLA and are called from fitabn().

}

\author{Fraser Ian Lewis}


\keyword{internal}
