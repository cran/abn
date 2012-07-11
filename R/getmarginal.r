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

getmarginal <- function(data.df,dag.m,whichnode=NULL, whichvar="constant", 
                        hyper.params=list("mean"=c(0),"sd"=c(sqrt(1000)),"shape"=c(0.001),"scale"=c(1/0.001)),
                        verbose=FALSE,post.x=seq(-0.1,0.1,len=100),
                        max.iters=100,epsabs=1e-7,std=TRUE) {
    
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
    data.df.cat<-data.df[,get.factors];
    obsdata.cat<-makeintofactors(data.df.cat);## need to make into integers
    numVarLevels<-as.integer(apply(obsdata.cat,2,max));## get vector of number of levels in each variables
    if( length(which(numVarLevels!=2)) ){#cat("Error in ",names(data.df)[which(min(numVarLevels)==1)],"\n");
                             stop("all categorical variables must have exactly TWO categories\n - create additional binary variables for multinomial data\n");}
    }   
    
    if(length(get.gaussian)>0){
                               obsdata.cts<-as.data.frame(data.df[,get.gaussian]);names(obsdata.cts)<-names(data.df)[get.gaussian];### NOTE as.data.frame() is NEW
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
                                                         
    tmp<-makeintobinarymatrix(dag.m); ## cohersion to ints and simple sanity checks
    dag<-tmp$dag;
    maxparents<-max(tmp$maxparents,1);## if a fully independent DAG then 1 is a dummy value required for memory allocation
    maxparents<-as.integer(maxparents);

    if(!is.logical(verbose)){stop("verbose must be either TRUE or FALSE");}
    if(verbose){loc.verbose<-1;} else {loc.verbose<-0;}
    loc.verbose<-as.integer(loc.verbose);

    loc.whichnode<-as.integer(which(rownames(dag.m)==whichnode)); ## which response variable e.g. which row in dag.m matrix
    if(length(loc.whichnode)==0){stop("unknown node name");}
    loc.whichgaus<- 0; # need index of gaussian variable e.g. is it the first or second gaussian variable etc
    if(var.types[loc.whichnode]==0){#have gaussian node
                      tmp<-var.types[1:loc.whichnode];
                      loc.whichgaus<-length(tmp[tmp==0]); if(loc.whichgaus<1){stop("error");}
                      }

    if(var.types[loc.whichnode]!=0 #not a gaussian node 
                                  && whichvar=="precision"){stop("precision can only be used with Gaussian nodes!");} 
                             
    ## note whichvar is the index from the left only counting cols with 1 e.g. model covariates
    if(whichvar=="constant" || whichvar=="precision"){
                                              if(whichvar=="constant"){loc.whichvariable<-as.integer(1);
                                              } else {tmp<-dag.m[loc.whichnode,];tmp<-tmp[tmp==1];
                                                      loc.whichvariable<-as.integer(length(tmp)+2);}#+2 as +1 for constant then another +1 for sd
    } else {tmp<-dag.m[loc.whichnode,];tmp<-tmp[tmp==1];got<-which(names(tmp)==whichvar)+1;#+1 since always have the constant term first
    loc.whichvariable<-as.integer(got);} ## which covariate e.g. which col in dag.m matrix
    if(length(loc.whichvariable)==0){stop("unknown variable name");}
    #cat("numnode=",loc.whichnode," numvar=",loc.whichvariable,"\n");
    #cat("whichgaus=",loc.whichgaus,"\n");
    
    loc.matrix.posterior<-cbind(as.real(post.x),as.real(rep(-1,length(post.x))));## create matrix with first col x, second col with be f(x)
    loc.numvariates<-dim(loc.matrix.posterior)[1];# number of times to calculate posterior
    
    ## coerce list into two vectors, one for means and one for **standard deviation** 
    prior.mean<-hyper.params$mean;
    if(length(prior.mean)!=dim(data.df)[2]+1){prior.mean<-as.double(rep(0.0,dim(data.df)[2]+1));}      
    #if(length(prior.mean)!=dim(data.df)[2]+1){stop("need prior mean for each variable in data.frame");}
    
    prior.sd  <-(hyper.params$sd);
    if(length(prior.sd)!=dim(data.df)[2]+1){prior.sd<-as.double(rep(sqrt(1000.0),dim(data.df)[2]+1));}      
    #if(length(prior.sd)!=dim(data.df)[2]+1){stop("need prior mean for each variable in data.frame");}
    
    prior.gamma.shape<-hyper.params$shape;
    if(length(get.gaussian)>0 && length(prior.gamma.shape)!=length(get.gaussian)){prior.gamma.shape<-as.double(rep(0.001,length(get.gaussian)));}  
    prior.gamma.scale<-hyper.params$scale;
    if(length(get.gaussian)>0 && length(prior.gamma.scale)!=length(get.gaussian)){prior.gamma.scale<-as.double(rep(1/0.001,length(get.gaussian)));}
    
    ## new part 31-jan-2012. To fix bug if only passed a discrete model R croaks about REAL() applied to non-numbers if prior.gamma.shape and prior.gamma.scale 
    ##        are not set to something. Use a dummy of zero.
    if(length(get.gaussian)==0){prior.gamma.shape<-0.0;prior.gamma.scale<-0.0;}

    max.iters<-as.integer(max.iters);
    epsabs<-as.double(epsabs);
    ## call to C
        res <- .Call("getmarginals_additive",obsdata,dag,prior.mean,prior.sd,prior.gamma.shape,prior.gamma.scale,maxparents,loc.verbose,var.types, 
                                       loc.matrix.posterior, loc.numvariates, loc.whichnode-1, loc.whichvariable-1,loc.whichgaus-1,max.iters,epsabs

               ,PACKAGE="abn" ## uncomment to load as package not shlib
              )
    res<-matrix(unlist(res),ncol=2);colnames(res)<-c("x","f");
    return(res);
             
}
   

