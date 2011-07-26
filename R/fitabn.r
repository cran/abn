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

fitabn <- function(data.df,dag.m,hyper.params=list("mean"=c(0),"var"=c(1000)), verbose=FALSE) {
    
    obsdata<-makeintofactors(data.df);## need to make into integers
    tmp<-makeintobinarymatrix(dag.m); ## cohersion to ints and simple sanity checks
    dag<-tmp$dag;
    maxparents<-max(tmp$maxparents,1);## if a fully independent DAG then 1 is a dummy value required for memory allocation
    maxparents<-as.integer(maxparents); 
    
    numVarLevels<-as.integer(apply(obsdata,2,max));## get vector of number of levels in each variables
    if( length(which(numVarLevels!=2)) ){#cat("Error in ",names(data.df)[which(min(numVarLevels)==1)],"\n");
                             stop("sorry, currently all variables in additive models must have exactly TWO categories\n");}
       

    if(!is.logical(verbose)){stop("verbose must be either TRUE or FALSE");}
    if(verbose){loc.verbose<-1;} else {loc.verbose<-0;}
    loc.verbose<-as.integer(loc.verbose);

    ## create a list containing the names of all the nodes and the levels in each node
    loc.labels<-list(names=names(obsdata));
    for(i in 1:(dim(obsdata)[2])){loc.labels[[i+1]]<-levels(as.factor(obsdata[,i]));} #use data.df since obsdata HAS NO LEVELS set
    names(loc.labels)[2:length(loc.labels)]<-names(obsdata);
    
    ## coerce list into two vectors, one for means and one for **standard deviation** 
    prior.mean<-hyper.params$mean;
    if(length(prior.mean)!=dim(data.df)[2]+1){prior.mean<-as.double(rep(0.0,dim(data.df)[2]+1));}      
    if(length(prior.mean)!=dim(data.df)[2]+1){stop("need prior mean for each variable in data.frame");}
    
    prior.sd  <-sqrt(hyper.params$var);
    if(length(prior.sd)!=dim(data.df)[2]+1){prior.sd<-as.double(rep(sqrt(1000.0),dim(data.df)[2]+1));}      
    if(length(prior.sd)!=dim(data.df)[2]+1){stop("need prior mean for each variable in data.frame");}
    #cat("Remember to constrain C code to a single node in laplace_score...()\n");
    ## call to C
    res <- .Call("fitnetwork_additive",obsdata,dag,prior.mean,prior.sd,maxparents,numVarLevels,loc.labels,loc.verbose
              ,PACKAGE="abn" ## uncomment to load as package not shlib
              )
    return(unlist(res));
}
   

