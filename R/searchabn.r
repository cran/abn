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

searchabn <- function(data.df,banned.m,hyper.params=list("mean"=c(0),"var"=c(1000)), max.parents=NULL,init.permuts=0) {
    
    obsdata<-makeintofactors(data.df);## need to make into integers
    tmp<-makeintobinarymatrix(banned.m); ## cohersion to ints and simple sanity checks
    dag<-tmp$dag;
    maxparents<-max(max.parents,1);## if a fully independent DAG then 1 is a dummy value required for memory allocation
    maxparents<-as.integer(maxparents); 
    
    numVarLevels<-as.integer(apply(obsdata,2,max));## get vector of number of levels in each variables
    if( length(which(numVarLevels!=2)) ){#cat("Error in ",names(data.df)[which(min(numVarLevels)==1)],"\n");
                             stop("sorry, currently all variables in additive models must have exactly TWO categories\n");}  

    ## create a list containing the names of all the nodes and the levels in each node
    loc.labels<-list(names=names(obsdata));
    for(i in 1:(dim(obsdata)[2])){loc.labels[[i+1]]<-levels(as.factor(obsdata[,i]));} #use data.df since obsdata HAS NO LEVELS set
    names(loc.labels)[2:length(loc.labels)]<-names(obsdata);
    
    nopermuts<-as.integer(init.permuts);
    ## we use a random network to start search heuristic and to create this we need a random shuffle - done in C using GSL previously
    maxnumlinks<-dim(data.df)[2]*(dim(data.df)[2]-1);## upper bound number of links ignoring cycles
    shuffle<-sample(0:(maxnumlinks-1));shuffle<-as.integer(shuffle);##C indexes start at 0!

    ## coerce list into two vectors, one for means and one for **standard deviation** 
    prior.mean<-as.double(hyper.params$mean);      
    if(length(hyper.params$mean)!=dim(data.df)[2]+1){stop("need prior mean for each variable in data.frame");}
    prior.sd  <-as.double(sqrt(hyper.params$var));
    if(length(hyper.params$var)!=dim(data.df)[2]+1){stop("need prior mean for each variable in data.frame");}
    ## call to C
    res <- .Call("searchfornetwork_additive",obsdata,dag,maxparents,prior.mean,prior.sd,numVarLevels,nopermuts,shuffle,loc.labels
               ,PACKAGE="abn" ## uncomment to load as package not shlib
              )
    for(i in 2:length(res)){colnames(res[[i]])<-names(obsdata);rownames(res[[i]])<-names(obsdata);}
    names(res)<-c("scores",paste("iter",1:(length(res)-1),sep=""));
    return((res));

}
   

