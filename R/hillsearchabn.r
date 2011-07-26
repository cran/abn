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

hillsearchabn <- function(data.df,banned.m,hyper.params=list("mean"=c(0),"var"=c(1000)), max.parents=NULL,init.permuts=0,num.searches=1) {
    
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
    nosearches<-as.integer(num.searches);

    ## we use a random network to start search heuristic and to create this we need a random shuffle - done in C using GSL previously
    maxnumlinks<-dim(data.df)[2]*(dim(data.df)[2]-1);## upper bound number of links ignoring cycles
    ## unlike searchbn() we pass all the random shuffles needed e.g. num.searches X 
    store<-matrix(rep(0:(maxnumlinks-1),num.searches),nrow=maxnumlinks);##use matrix each COL will be a shuffle
    for(i in 1:num.searches){store[,i]<-sample(0:(maxnumlinks-1));} ## this avoid slow store<-c(store,) etc
    shuffle<-as.vector(store);## turn into one long array
    shuffle<-as.integer(shuffle);## e.g. shuffle[1:maxlinks] is first shuffle, shuffle[maxlinks+1:2*maxlinks] is second etc

    ## coerce list into two vectors, one for means and one for **standard deviation** 
    prior.mean<-as.double(hyper.params$mean);      
    if(length(hyper.params$mean)!=dim(data.df)[2]+1){stop("need prior mean for each variable in data.frame");}
    prior.sd  <-as.double(sqrt(hyper.params$var));
    if(length(hyper.params$var)!=dim(data.df)[2]+1){stop("need prior mean for each variable in data.frame");}
    ## call to C
    res <- .Call("hillsearchfornetwork_additive",obsdata,dag,maxparents,prior.mean,prior.sd,numVarLevels,nopermuts,shuffle,loc.labels, nosearches
               ,PACKAGE="abn" ## uncomment to load as package not shlib
              )
    ## results format is res[[1]] is a vector of network scores in order: init score1, final score1, init score2, final score2,....etc
    ## with res[[2]] = init network 1, res[[3]] final network 1, res[[4]] init network2, res[[5]] init network2,... etc
    for(i in 2:length(res)){colnames(res[[i]])<-names(obsdata);rownames(res[[i]])<-names(obsdata);}
    ## now re-organise res into something easier to analyse - a list of three lists
    ## list 1 - vector of scores and list of initial matrices, list 2 -vector of score and list of final matrices
    scores<-res[[1]];
     init.indexes<-seq(1,2*num.searches,by=2);
     fin.indexes<-seq(2,2*num.searches,by=2);
                     
    init.scores<-scores[init.indexes];## scores from initial networks
    fin.scores<-scores[fin.indexes];
    init.mat<-res[init.indexes+1];##offset for score vector
    fin.mat<-res[fin.indexes+1];##offset for score vector
    rm(res);

    return(list(init.score=init.scores,final.score=fin.scores,init.dag=init.mat,final.dag=fin.mat));

}
   

