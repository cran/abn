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

findmostprobablebn <- function(local.scores, data.df=NULL) {
    
   #if(local.scores$restricted.parents==1 && feature != "custom"){stop("Must use custom featuredefn when a parent limit has been imposed (inc. for denominator)");}

    loc.numnodes<-as.integer(dim(local.scores$parents)[2]);
    loc.maxparents<-max(apply(local.scores$parents,1,sum));#maximum number of parents in any node
    local.scores$node<-as.integer(local.scores$node-1);#since C indexes from 0

    #if(is.null(offset)){offset<-exp(-mean(local.scores$nodescore));} #use mean as seems to work fine n.b. =1 does not work for larger sample sizes
    offset<-1.0;
    offset<-as.double(offset);
    #offset<-exp(-mean(local.scores$nodescore));##to help avoid underflow - mult.factor on parent prior
    ## call to C
    res.prob <- .Call("getposterior_features_max",local.scores,loc.numnodes,loc.maxparents,offset
               ,PACKAGE="abn" ## uncomment to load as package not shlib
              );
    junk<-gc(FALSE);## some garbage collection 
  
   if(TRUE){  res<-list();
     tmpmat<-matrix(rep(-1,loc.numnodes*loc.numnodes),ncol=loc.numnodes);
     #colnames(tmpmat)<-names(data.df);rownames(tmpmat)<-names(data.df);
     tmpvec<-rep(-1,loc.numnodes);
     res[[1]]<-tmpvec;res[[2]]<-tmpmat;
     for(i in 1:length(res.prob)){res[[1]][i]<-res.prob[[i]][1]+1;
                                  res[[2]][i,]<-res.prob[[i]][-1];}
     names(res[[1]])<-1:length(res[[1]]);
     srt<-as.numeric(names(sort(res[[1]])));
     res[[2]]<-res[[2]][srt,];
     colnames(res[[2]])<-names(data.df);rownames(res[[2]])<-names(data.df);
    
    return(res[[2]]);
}
#return(res.prob);
    #loc.child<- -1; loc.child<-as.integer(loc.child);
    #res.denom <- .Call("getposterior_features",local.scores,loc.numnodes,-1,loc.parent,loc.maxparents
               #,PACKAGE="abn" ## uncomment to load as package not shlib
     #         );
    #junk<-gc(FALSE);## some garbage collection 
    #res[[2]]<-matrix(res[[2]],ncol=dim(data.df)[2],byrow=TRUE);## need to tranform long vector into matrix of parent combinations 
    #names(res)<-c("node","parents","nodescore");
    
    ## res[[1]] is a vector of node indexes (from 1 not zero)
    ## res[[2]] is a matrix where each row is a parent combination for the node in res[[1]]
    ## res[[3]] is the score for the parent combination in res[[2]]
}
   

