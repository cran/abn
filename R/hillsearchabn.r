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

hillsearchabn <- function(data.df,banned.m,retain.m,start.m,hyper.params=list("mean"=c(0),"sd"=c(sqrt(1000)),"shape"=c(0.001),"scale"=c(1/0.001)),
                          max.parents=NULL,init.permuts=0,num.searches=1,db.size=10000,localdb=TRUE,timing=TRUE,max.iters=100,epsabs=1e-7,error.verbose=FALSE,enforce.db.size=TRUE) {
    
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
                               obsdata.cts<-data.df[,get.gaussian];
                               for(i in 1:dim(obsdata.cts)[2]){obsdata.cts[,i]<-as.double(obsdata.cts[,i]);}
    }
    
    #now put all the data back into a single data frame
    if(length(get.factors)>0 && length(get.gaussian)==0){obsdata<-data.frame(obsdata.cat);} #categorical only
    if(length(get.factors)==0 && length(get.gaussian)>0){obsdata<-data.frame(obsdata.cts);} #gaussian only
    if(length(get.factors)>0 && length(get.gaussian)>0){ obsdata<-data.frame(obsdata.cat,obsdata.cts); #categorical and gaussian
                                                         obsdata<-obsdata[,names(data.df)];           } #get into original colorder
                                                      
    for(i in 1:dim(obsdata)[2]){obsdata[,i]<-as.double(obsdata[,i]);}#coerce ALL cols to double equivalents
                                                         
    dag<-makeintobinarymatrix(banned.m,ban=TRUE)$dag; ## cohersion to ints and simple sanity checks
    dag.retain<-makeintobinarymatrix(retain.m,ban=FALSE)$dag;                                                                                                             
    if(!is.list(start.m)){stop("start.m must be a list");}
    if(length(start.m)!=num.searches){stop("start.m must same length as num.searches");}
    dag.start<-list();
    for(i in 1:num.searches){dag.start[[i]]<-makeintobinarymatrix(start.m[[i]],ban=FALSE)$dag;} #get a list of validated matrices
    db.size<-as.integer(db.size);
    maxparents<-max(max.parents,1);## if a fully independent DAG then 1 is a dummy value required for memory allocation
    maxparents<-as.integer(maxparents); 
    if(max(apply(retain.m,1,sum))>maxparents){stop("retain.m is inconsistent with max.parents");}
    #numVarLevels<-as.integer(apply(obsdata,2,max));## get vector of number of levels in each variables
    #if( length(which(numVarLevels!=2)) ){#cat("Error in ",names(data.df)[which(min(numVarLevels)==1)],"\n");
    #                         stop("sorry, currently all variables in additive models must have exactly TWO categories\n");}  

    ## create a list containing the names of all the nodes and the levels in each node
    #loc.labels<-list(names=names(obsdata));
    #for(i in 1:(dim(obsdata)[2])){loc.labels[[i+1]]<-levels(as.factor(obsdata[,i]));} #use data.df since obsdata HAS NO LEVELS set
    #names(loc.labels)[2:length(loc.labels)]<-names(obsdata);
    if(!is.logical(timing)){stop("timing must be either TRUE or FALSE");}
    if(timing){loc.timing<-1;## turn on timing output
    } else {loc.timing<-0;}
    if(!is.logical(localdb)){stop("localdb must be either TRUE or FALSE");}
    if(localdb){loc.localdb<-1;## turn on timing output
    } else {loc.localdb<-0;}
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
    if(length(prior.mean)!=dim(data.df)[2]+1){prior.mean<-as.double(rep(0.0,dim(data.df)[2]+1));
                                              warning("using prior mean of 0.0 for each variable");}
    prior.sd  <-as.double(hyper.params$sd);
    if(length(prior.sd)!=dim(data.df)[2]+1){prior.sd<-as.double(rep(sqrt(1000.0),dim(data.df)[2]+1));
                                            warning("using prior sd of sqrt(1000) for each variable");}
    
    if(is.null(hyper.params$shape)){hyper.params$shape=0.001;}
    prior.gamma.shape<-hyper.params$shape;
    if(length(get.gaussian)>0 && length(prior.gamma.shape)!=length(get.gaussian)){prior.gamma.shape<-as.double(rep(0.001,length(get.gaussian)));
               warning("using 0.001 as prior shape for each gaussian variable in data.frame"); }
    if(is.null(hyper.params$scale)){hyper.params$scale=1/0.001;}  
    prior.gamma.scale<-hyper.params$scale;
    if(length(get.gaussian)>0 && length(prior.gamma.scale)!=length(get.gaussian)){prior.gamma.scale<-as.double(rep(1/0.001,length(get.gaussian)));
               warning("using 1/0.001 as prior scale for each gaussian variable in data.frame");}
    max.iters<-as.integer(max.iters);
    epsabs<-as.double(epsabs);
    
    if(!is.logical(error.verbose)){stop("erorr.verbose must be either TRUE or FALSE");}
    if(error.verbose){loc.error.verbose<-1;} else {loc.error.verbose<-0;}
    loc.error.verbose<-as.integer(loc.error.verbose);
    
    if(!is.logical(enforce.db.size)){stop("enforce.db.size must be either TRUE or FALSE");}
    if(enforce.db.size){loc.enforce.db.size<-1;} else {loc.enforce.db.size<-0;}
    loc.enforce.db.size<-as.integer(loc.enforce.db.size);
    
    ## call to C
    res <- .Call("hillsearchfornetwork_additive",obsdata,dag,maxparents,prior.mean,prior.sd,prior.gamma.shape,prior.gamma.scale,var.types,
                                                nopermuts,shuffle, nosearches,dag.retain,dag.start,db.size,loc.localdb,loc.timing,max.iters,epsabs,
                                                loc.error.verbose, loc.enforce.db.size
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
   

