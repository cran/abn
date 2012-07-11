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

allnodesabn <- function(data.df,hyper.params=list("mean"=c(0),"sd"=c(sqrt(1000)),"shape"=c(0.001),"scale"=c(1/0.001)),max.iters=100,epsabs=1e-7, 
                   verbose=FALSE, error.verbose=FALSE, max.parents=NULL, all.nodes=TRUE, which.nodes=NULL,std=TRUE) {
    
    #need to find categorical variables and continuous - Gaussian
    get.factors<-NULL; 
    get.gaussian<-NULL;
    var.types<-rep(-1,dim(data.df)[2]);
    for(i in 1:dim(data.df)[2]){#find which cols are which variable type
                                if(is.factor(data.df[,i])){get.factors<-c(get.factors,i);
                                   var.types[i]<-1;#got categorical variable
                                } else {get.gaussian<-c(get.gaussian,i);
                                   var.types[i]<-0;# got numeric variable
                                   }
    }
    var.types<-as.integer(var.types);#1-cat, 0-numeric
    
    if(length(get.factors)>0){#have some categorical variables so some coersion to factors then integers needed
    data.df.cat<-data.frame(data.df[,get.factors]);names(data.df.cat)<-names(data.df)[get.factors];
    obsdata.cat<-makeintofactors(data.df.cat);## need to make into integers
    numVarLevels<-as.integer(apply(obsdata.cat,2,max));## get vector of number of levels in each variables
    if( length(which(numVarLevels!=2)) ){#cat("Error in ",names(data.df)[which(min(numVarLevels)==1)],"\n");
                             stop("all categorical variables must have exactly TWO categories\n - create additional binary variables for multinomial data\n");}
    }   
    
    if(length(get.gaussian)>0){
                               obsdata.cts<-data.frame(data.df[,get.gaussian]);names(obsdata.cts)<-names(data.df)[get.gaussian];
                               for(i in 1:dim(obsdata.cts)[2]){if(std){# std. gaus vars to mean zero and sd=1
                                                               obsdata.cts[,i]<- (obsdata.cts[,i]-mean(obsdata.cts[,i]))/sd(obsdata.cts[,i]);}
                                                               obsdata.cts[,i]<-as.double(obsdata.cts[,i]);}
    }
    
    #now put all the data back into a single data frame
    if(length(get.factors)>0 && length(get.gaussian)==0){obsdata<-data.frame(obsdata.cat);} #categorical only
    if(length(get.factors)==0 && length(get.gaussian)>0){obsdata<-data.frame(obsdata.cts);} #gaussian only
    if(length(get.factors)>0 && length(get.gaussian)>0){ obsdata<-data.frame(obsdata.cat,obsdata.cts); #categorical and gaussian
                                                         obsdata<-obsdata[,names(data.df)];           } #get into original colorder
                                                      
    for(i in 1:dim(obsdata)[2]){obsdata[,i]<-as.double(obsdata[,i]);}#coerce ALL cols to double equivalents
                                                         
    maxparents<-max(max.parents,1);## if a fully independent DAG then 1 is a dummy value required for memory allocation
    if(maxparents>=dim(data.df)[2]){stop("max.parents must be strictly less than the number of nodes in the network!");}
    maxparents<-as.integer(maxparents); 
     
    if(!is.logical(verbose)){stop("verbose must be either TRUE or FALSE");}
    if(verbose){loc.verbose<-1;} else {loc.verbose<-0;}
    loc.verbose<-as.integer(loc.verbose);

    if(!is.logical(error.verbose)){stop("error.verbose must be either TRUE or FALSE");}
    if(error.verbose){loc.error.verbose<-1;} else {loc.error.verbose<-0;}
    loc.error.verbose<-as.integer(loc.error.verbose);

    ## coerce list into two vectors, one for means and one for **standard deviation** 
    prior.mean<-hyper.params$mean;
    if(length(prior.mean)!=dim(data.df)[2]+1){prior.mean<-as.double(rep(0.0,dim(data.df)[2]+1));}      
    #if(length(prior.mean)!=dim(data.df)[2]+1){stop("need prior mean for each variable in data.frame");}
    
    prior.sd  <-hyper.params$sd;
    if(length(prior.sd)!=dim(data.df)[2]+1){prior.sd<-as.double(rep(sqrt(1000.0),dim(data.df)[2]+1));}      
    #if(length(prior.sd)!=dim(data.df)[2]+1){stop("need prior mean for each variable in data.frame");}
    
    #if(is.null(hyper.params$shape)){hyper.params$shape=0.001;}
    prior.gamma.shape<-hyper.params$shape;
    if(length(get.gaussian)>0 && length(prior.gamma.shape)!=length(get.gaussian)){prior.gamma.shape<-as.double(rep(0.001,length(get.gaussian)));}
    #if(is.null(hyper.params$scale)){hyper.params$scale=1/0.001;}  
    prior.gamma.scale<-hyper.params$scale;
    if(length(get.gaussian)>0 && length(prior.gamma.scale)!=length(get.gaussian)){prior.gamma.scale<-as.double(rep(1/0.001,length(get.gaussian)));}
    
    max.iters<-as.integer(max.iters);
    epsabs<-as.double(epsabs);

    ## which.nodes to consider - to allow computation to be split over different cpus if need be
    if(!is.logical(all.nodes)){stop("all.nodes must be either TRUE or FALSE");} 
    if(all.nodes){which.nodes<-1:dim(data.df)[2];## all nodes
    } else { if(min(which.nodes)<1 || max(which.nodes)> dim(data.df)[2]){stop("which.nodes is invalid!");}}
    which.nodes<-as.integer(which.nodes); 
    
    #cat("Remember to constrain C code to a single node in laplace_score...()\n");
    ## call to C
    res <- .Call("allnodesbn_additive",obsdata,prior.mean,prior.sd,prior.gamma.shape,prior.gamma.scale,maxparents,loc.verbose,var.types,max.iters,epsabs,
                                        loc.error.verbose,which.nodes
              ,PACKAGE="abn" ## uncomment to load as package not shlib
              )  
    res[[2]]<-matrix(res[[2]],ncol=dim(data.df)[2],byrow=TRUE);## need to tranform long vector into matrix of parent combinations 
    names(res)<-c("node","parents","nodescore","restricted.parents");
    return(res);        
   
}
   

