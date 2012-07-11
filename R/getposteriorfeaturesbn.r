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

getposteriorfeaturesbn <- function(local.scores,feature="arc",child=NULL, parent=NULL, featuredefn=NULL, offset=1) {
    
   #if(local.scores$restricted.parents==1 && feature != "custom"){stop("Must use custom featuredefn when a parent limit has been imposed (inc. for denominator)");}

   if(feature=="arc" && !is.null(child) && !is.null(parent) ){loc.child<-as.integer(child-1);#since C indexes from 0
                                                               loc.parent<-as.integer(parent-1);#since C indexes from 0
                                                               featuredefn<-NULL; #avoid errors in case this is given 
    } else { if(feature=="custom"){ #user-supplied definition
                                    if(length(featuredefn)!=dim(local.scores$parents)[1]){stop("incorrect1 format for featuredefn");}
                                    featuredefn<-as.integer(featuredefn);
                                    if(length(which(featuredefn==1))+
                                       length(which(featuredefn==0))!=length(featuredefn) ){stop("incorrect2 format for featuredefn"); }
                                    loc.child<- NULL;
                                    loc.parent<-NULL;   
   
            }  else {if(feature=="all"){loc.child<- -1;## a simply flag to say that all features are included
                                                       ## this is used to get denominator
                                        loc.child<-as.integer(loc.child);
                                        loc.parent<-as.integer(parent-1);# A DUMMY - note used
                                        featuredefn<-NULL;#avoid errors in case this is given
                    
                      } else {stop("unknown type of feature");}}}
    loc.numnodes<-as.integer(dim(local.scores$parents)[2]);
    loc.maxparents<-max(apply(local.scores$parents,1,sum));#maximum number of parents in any node
    local.scores$node<-as.integer(local.scores$node-1);#since C indexes from 0
   # if(loc.maxparents<(loc.numnodes-1)){stop("maxparents must equal numnodes-1: code not yet completed for parent limit k < numnodes-1!");}

    #if(is.null(offset)){offset<-exp(-mean(local.scores$nodescore));} #use mean as seems to work fine n.b. =1 does not work for larger sample sizes
    offset<-as.double(offset);
    #offset<-exp(-mean(local.scores$nodescore));##to help avoid underflow - mult.factor on parent prior
    ## call to C
    res.prob <- .Call("getposterior_features",local.scores,loc.numnodes,loc.child,loc.parent,loc.maxparents,featuredefn,offset
               ,PACKAGE="abn" ## uncomment to load as package not shlib
              );
    junk<-gc(FALSE);## some garbage collection 
  
    return(res.prob[[1]][1]);

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
   

