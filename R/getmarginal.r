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

getmarginal <- function(data.df,dag.m,whichnode=NULL, whichvar="constant", hyper.params=list("mean"=c(0),"var"=c(1000)),
                        verbose=FALSE,post.x=seq(-0.1,0.1,len=100)) {
    
    if(length(rownames(dag.m)) && length(colnames(dag.m))){ } else stop("need row and col names set");
    if(length(which(colnames(dag.m)==whichnode))==0){stop("invalid whichnode name");}
    if(length(which(colnames(dag.m)==whichvar))==0 && whichvar!="constant"){stop("invalid whichvar name");}

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

    loc.whichnode<-as.integer(which(rownames(dag.m)==whichnode)); ## which response variable e.g. which row in dag.m matrix
    ## note whichvar is the index from the left only counting cols with 1 e.g. model covariates
    if(whichvar=="constant"){loc.whichvariable<-as.integer(1);
    } else {tmp<-dag.m[loc.whichnode,];tmp<-tmp[tmp==1];got<-which(names(tmp)==whichvar)+1;#+1 since always have the constant term first
    loc.whichvariable<-as.integer(got);} ## which covariate e.g. which col in dag.m matrix
    #cat("numnode=",loc.whichnode," numvar=",loc.whichvariable,"\n");
    
    loc.matrix.posterior<-cbind(as.real(post.x),as.real(rep(-1,length(post.x))));## create matrix with first col x, second col with be f(x)
    loc.numvariates<-dim(loc.matrix.posterior)[1];# number of times to calculate posterior
    ## create a list containing the names of all the nodes and the levels in each node
    loc.labels<-list(names=names(obsdata));
    for(i in 1:(dim(obsdata)[2])){loc.labels[[i+1]]<-levels(as.factor(obsdata[,i]));} #use data.df since obsdata HAS NO LEVELS set
    names(loc.labels)[2:length(loc.labels)]<-names(obsdata);
    
    ## coerce list into two vectors, one for means and one for **standard deviation** 
    prior.mean<-as.double(hyper.params$mean);      
    if(length(hyper.params$mean)!=dim(data.df)[2]+1){stop("need prior mean for each variable in data.frame");}
    prior.sd  <-as.double(sqrt(hyper.params$var));
    if(length(hyper.params$var)!=dim(data.df)[2]+1){stop("need prior mean for each variable in data.frame");}
    ## call to C
    res <- .Call("getmarginals_additive",obsdata,dag,prior.mean,prior.sd,maxparents,numVarLevels,loc.labels,loc.verbose, 
                                       loc.matrix.posterior, loc.numvariates, loc.whichnode-1, loc.whichvariable-1
               ,PACKAGE="abn" ## uncomment to load as package not shlib
              )
    res<-matrix(unlist(res),ncol=2);colnames(res)<-c("x","f");
    return(res);
}
   

