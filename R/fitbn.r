## fitbn.R --- 
## Author          : Fraser Lewis
## Created On      : Sun May 13:43 2010
## Last Modified By: Fraser Lewis
## Last Modified On: Sun May 13:43 2010
## Update Count    : 0
## Status          : Unknown, Use with caution!
###############################################################################
##
##    Copyright (C) 2002  Fraser Lewis
##
##    This program is free software; you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation; either version 2 of the License, or
##    (at your option) any later version.
##
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with this program; if not, write to the Free Software
##    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
######################################################################

fitbn <- function(data.df,dag.m,prior.obs.per.node=NULL,useK2=FALSE, verbose=FALSE) {
    
    obsdata<-makeintofactors(data.df);## need to make into integers
    tmp<-makeintobinarymatrix(dag.m); ## cohersion to ints and simple sanity checks
    dag<-tmp$dag;
    maxparents<-max(tmp$maxparents,1);## if a fully independent DAG then 1 is a dummy value required for memory allocation
    maxparents<-as.integer(maxparents); 
    if(!is.logical(useK2)){stop("useK2 must be either TRUE or FALSE");}
    if(useK2){loc.useK2<-1;## use K2
    } else {loc.useK2<-0;  ## use BDeu
            if(is.numeric(prior.obs.per.node) && prior.obs.per.node>0.0){
              } else {stop("invalid prior.obs.per.node - must be numeric and strictly positive");}
            }
    if(!is.logical(verbose)){stop("verbose must be either TRUE or FALSE");}
    if(verbose){loc.verbose<-1;} else {loc.verbose<-0;}
    loc.verbose<-as.integer(loc.verbose);
    numVarLevels<-as.integer(apply(obsdata,2,max));## get vector of number of levels in each variables
    if(min(numVarLevels)==1){#cat("Error in ",names(data.df)[which(min(numVarLevels)==1)],"\n");
                             stop("variables must have at least TWO categories - ",names(data.df)[which(numVarLevels==1)],"- does not\n");}
    prior.obs.per.node<-as.double(prior.obs.per.node);## coerce just in case
    ## create a list containing the names of all the nodes and the levels in each node
    loc.labels<-list(names=names(obsdata));
    for(i in 1:(dim(obsdata)[2])){loc.labels[[i+1]]<-levels(as.factor(obsdata[,i]));} #use data.df since obsdata HAS NOT LEVELS set
    names(loc.labels)[2:length(loc.labels)]<-names(obsdata);
    ## call to C
    res <- .Call("fitnetwork",obsdata,dag,loc.useK2,maxparents,prior.obs.per.node,numVarLevels,loc.labels,loc.verbose
               ,PACKAGE="abn" ## uncomment to load as package not shlib
              )
    return(unlist(res));
}
   

