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

searchbn <- function(data.df,banned.m,prior.obs.per.node=NULL,useK2=FALSE,max.parents=NULL,init.permuts=0) {
    
    obsdata<-makeintofactors(data.df);
    tmp<-makeintobinarymatrix(banned.m);
    dag<-tmp$dag;
    maxparents<-max(max.parents,1);## if a fully independent DAG then 1 is a dummy value required for memory allocation
    maxparents<-as.integer(maxparents); 
    if(!is.logical(useK2)){stop("useK2 must be either TRUE or FALSE");}
    if(useK2){loc.useK2<-1;## use K2
    } else {loc.useK2<-0;  ## use BDeu
            if(is.numeric(prior.obs.per.node) && prior.obs.per.node>0.0){
              } else {stop("invalid prior.obs.per.node - must be numeric and strictly positive");}
            }

    ## create a list containing the names of all the nodes and the levels in each node -NOT USED but keeps network_score() happy
    loc.labels<-list(names=names(obsdata));
    for(i in 1:(dim(obsdata)[2])){loc.labels[[i+1]]<-levels(as.factor(obsdata[,i]));} #use data.df since obsdata HAS NOT LEVELS set
    names(loc.labels)[2:length(loc.labels)]<-names(obsdata);

    numVarLevels<-as.integer(apply(obsdata,2,max));#get vector of number of levels in each variables
    prior.obs.per.node<-as.double(prior.obs.per.node);#coerce just in case
    nopermuts<-as.integer(init.permuts);
    ## we use a random network to start search heuristic and to create this we need a random shuffle - done in C using GSL previously
    maxnumlinks<-dim(data.df)[2]*(dim(data.df)[2]-1);## upper bound number of links ignoring cycles
    shuffle<-sample(0:(maxnumlinks-1));shuffle<-as.integer(shuffle);##C indexes start at 0!
    ## call to C
    res <- .Call("searchfornetwork",obsdata,dag,loc.useK2,maxparents,prior.obs.per.node,numVarLevels,nopermuts,shuffle, loc.labels
              ,PACKAGE="abn" ## uncomment to load as package not shlib
              )
    for(i in 2:length(res)){colnames(res[[i]])<-names(obsdata);rownames(res[[i]])<-names(obsdata);}
    names(res)<-c("scores",paste("iter",1:(length(res)-1),sep=""));
    return((res));
}
   

