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

allnodesbn <- function(data.df,prior.obs.per.node=NULL,useK2=FALSE,max.parents=NULL, all.nodes=TRUE, which.nodes=NULL, verbose=FALSE) {
    
    obsdata<-makeintofactors(data.df);## need to make into integers
    
    maxparents<-max(max.parents,1);## if a fully independent DAG then 1 is a dummy value required for memory allocation
    if(maxparents>=dim(data.df)[2]){stop("max.parents must be strictly less than the number of nodes in the network!");}
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

    ## which.nodes to consider - to allow computation to be split over different cpus if need be
    if(!is.logical(all.nodes)){stop("all.nodes must be either TRUE or FALSE");} 
    if(all.nodes){which.nodes<-1:dim(data.df)[2];## all nodes
    } else { if(min(which.nodes)<1 || max(which.nodes)> dim(data.df)[2]){stop("which.nodes is invalid!");}}
    which.nodes<-as.integer(which.nodes); 
    
    ## call to C
    res <- .Call("allnodesbn",obsdata,loc.useK2,maxparents,prior.obs.per.node,numVarLevels,loc.verbose,which.nodes
               ,PACKAGE="abn" ## uncomment to load as package not shlib
              )
    res[[2]]<-matrix(res[[2]],ncol=dim(data.df)[2],byrow=TRUE);## need to tranform long vector into matrix of parent combinations 
    names(res)<-c("node","parents","nodescore","restricted.parents");
    return(res);
    ## res[[1]] is a vector of node indexes (from 1 not zero)
    ## res[[2]] is a matrix where each row is a parent combination for the node in res[[1]]
    ## res[[3]] is the score for the parent combination in res[[2]]
}
   

