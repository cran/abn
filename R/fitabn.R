###############################################################################
## fitabn.R --- 
## Author          : Fraser Lewis
## Last modified   : 29/08/2012
###############################################################################

## fit a given DAG to data using the distribution given
fitabn <- function(dag.m=NULL, data.df=NULL, data.dists=NULL, group.var=NULL,cor.vars=NULL,create.graph=FALSE,compute.fixed=FALSE,
                   return.modes=TRUE,
                   mean=0, prec=0.001,
                   loggam.shape=1,loggam.inv.scale=5e-05,verbose=FALSE, centre=TRUE,
                   max.iters=100,epsabs=1e-7,error.verbose=FALSE,epsabs.inner=1e-6,max.iters.inner=100,
                   finite.step.size=1E-07,hessian.params=c(1E-04,1E-02),max.iters.hessian=10,
                   min.pdf=1E-04,n.grid=NULL,
                   marginal.node=NULL, marginal.param=NULL,variate.vec=NULL
                  ){

      use.inla<-FALSE;
      return.modes<-TRUE;
      ntrials<-NULL;
      exposure<-NULL;
      ## check grouping variables
      list.group.var<-check.groups(group.var,data.df,cor.vars,use.inla);## returns ammended data.df and suitable variables
      data.df<-list.group.var$data.df;## this has removes the grouping variable from data.df
      grouped.vars<-list.group.var$grouped.vars;## int vect of variables to be treated as grouped
      group.ids<-list.group.var$group.ids;## int vector of group membership ids
      #data.df.inla<-list.group.var$data.df.inla;

      #print(data.df);
      ## run a series of checks on the data and distributions passed
      mylist<-check.data(data.df,data.dists,ntrials,exposure,use.inla,group.var);## return a list with entries bin, gaus, pois, ntrials and exposure

      ## run a series of checks on the DAG passed
      check.dag(dag.m,data.df=data.df,is.ban.matrix=FALSE,use.inla,group.var);

      ## coerce binary factors to become 0/1 integers - the 0 is based on the first entry in levels()
      if(!is.null(mylist$bin)){## have at least one binary variable
        for(i in mylist$bin){data.df[,i]<-as.numeric(data.df[,i])-1;}
      }

      ## standardize gaussian variables to zero mean and sd=1
      if(centre && !is.null(mylist$gaus)){## have at least one gaussian variable
        for(i in mylist$gaus){data.df[,i]<-(data.df[,i]-mean(data.df[,i]))/sd(data.df[,i]);}
      }

      ## check for type of grouped variable
      if(!is.null(group.var)){
      for(i in grouped.vars){## for each variable to be treated as grouped 
        if(data.dists[[i+1]]!="binomial" && !use.inla){## +1 is because grouped.vars is indexes from 0 for sending to C
                           stop("currently grouped variables must be binary only when use.inla=FALSE");}
      }}     
      ## Two computational options C or via calls to INLA

      res.list<-list();

      #########################################################
      ## use C
      #########################################################
      for(i in 1:dim(data.df)[2]){data.df[,i]<-as.double(data.df[,i]);}##coerce ALL cols to double equivalents
      var.types<-get.var.types(data.dists); ## get distributions in terms of a code
      ## C code only does limited number of dists relative to INLA - check that here
      if(length(which(var.types%in%c(4)))!=0){stop("This model cannot be fitted using fitabn()");}
      max.parents<-max(apply(dag.m,1,sum));
      #stop("");
      res <- .Call("fitabn",data.df,as.integer(dag.m),as.integer(dim(data.df)[2]),
                            as.integer(var.types),as.integer(max.parents),
                            as.double(mean),as.double(1/sqrt(prec)),as.double(loggam.shape),as.double(1/loggam.inv.scale),
                                         as.integer(max.iters),as.double(epsabs),
                                         as.integer(verbose),as.integer(error.verbose),
                             grouped.vars,## int.vector of variables which are mixed model nodes
                             group.ids,
                             as.double(epsabs.inner),
                             as.integer(max.iters.inner),
                             as.double(finite.step.size),
                             as.double(hessian.params),
                             as.integer(max.iters.hessian)
              ,PACKAGE="abn" ## uncomment to load as package not shlib
              )   
      ## now unroll results from C into a list
      for(child in 1:dim(dag.m)[1]){ 
                   child.name<-colnames(dag.m)[child];
                   res.list[[child.name]]<-res[[child]][1];}
      res.list[["mlik"]]<-sum(unlist(res.list));## add in total mlik for DAG
      
      res.list[["error.code"]]<-NULL;res.list[["hessian.accuracy"]]<-NULL;
      for(child in 1:dim(dag.m)[1]){
         res.list[["error.code"]]<-c(res.list[["error.code"]],res[[child]][2]);
         res.list[["hessian.accuracy"]]<-c(res.list[["hessian.accuracy"]],res[[child]][3]);}
         names(res.list[["error.code"]])<-names(res.list[["hessian.accuracy"]])<-colnames(dag.m);
      res.list[["error.code.desc"]]<-as.character(res.list[["error.code"]]);
      res.list[["error.code.desc"]]<-ifelse(res.list[["error.code.desc"]]==0,"success",res.list[["error.code"]]);
      res.list[["error.code.desc"]]<-ifelse(res.list[["error.code.desc"]]==1,"warning: mode results may be unreliable (optimiser terminated unusually)",res.list[["error.code.desc"]]);
      res.list[["error.code.desc"]]<-ifelse(res.list[["error.code.desc"]]==2,"error - logscore is NA - model could not be fitted",res.list[["error.code.desc"]]);
      #res.list[["error.code.desc"]]<-ifelse(res.list[["error.code.desc"]]==3,"warning: fin.diff adaptive step size failed (hessian NaN or Inf) - h_guess used",res.list[["error.code.desc"]]);
      res.list[["error.code.desc"]]<-ifelse(res.list[["error.code.desc"]]==4,"warning: fin.diff hessian estimation terminated unusually ",res.list[["error.code.desc"]]);

      #res.list[["use.inla"]]<-FALSE;
      
      if(return.modes){### Want the modes returned as part of the results
                      mymodes<-res;## make a copy
                      names(mymodes)<-colnames(dag.m);
                      for(i in 1:dim(dag.m)[1]){## for each node
                        mymodes[[i]]<-mymodes[[i]][-c(1:3)];## remove mlik - this is first entry, and error code and hessian accuracy
                        mymodes[[i]]<-mymodes[[i]][which(mymodes[[i]]!=.Machine$double.xmax)];## this discards all "empty" parameters 
                        nom<-colnames(dag.m)[which(dag.m[i,]==1)];
                        if(var.types[i]=="1"){nom<-c("(Intercept)",nom,"group.precision");}## binomial : just some naming for use later
                        if(var.types[i]=="2"){nom<-c("(Intercept)",nom,"precision","group.precision");} ## gaus
                        if(var.types[i]=="3"){nom<-c("(Intercept)",nom,"group.precision");} ## pois}
                        mynom<-NULL;
                        for(j in 1:length(mymodes[[i]])){mynom<-c(mynom,paste(colnames(dag.m)[i],nom[j],sep="|"));}
                        names(mymodes[[i]])<-mynom;                       
                      }
      
      res.list[["modes"]]<-mymodes;
      
      }
      #print(res); ## uncomment to print modes
      #####
      ##### Additional part only run if user wants marginal distributions
      #####
      if(compute.fixed){ res.list<-getmarginals(res.list, ## rest of arguments as for call to C fitabn
                                                data.df,dag.m,var.types,max.parents,
                                                mean,prec,loggam.shape,loggam.inv.scale,
                                                max.iters,epsabs,verbose,error.verbose,
                                                grouped.vars,## int.vector of variables which are mixed model nodes
                                                group.ids,
                                                epsabs.inner,max.iters.inner,finite.step.size,hessian.params,max.iters.hessian,
                                                min.pdf,marginal.node,marginal.param,variate.vec,n.grid);

        } ## end of compute.fixed 

      #########################################################
      ## Rgraph/graphviz part
      #########################################################
      if(create.graph){
      if(!require(Rgraphviz)){stop("library Rgraphviz is not available!\nRgraphviz is available as part of the bioconductor project - see http://www.bioconductor.org/install\nRgraphviz is required is create.graph=TRUE");}
      
      mygraph<-new("graphAM",adjMat=t(dag.m),edgemode="directed");
      res.list[["graph"]]<-mygraph;

      }

return(res.list);

}

#################################################################################################
## function for computing marginal posterior densities using C and is called as part of the 
## above fitabn() function. Only to be called internally.
#################################################################################################

getmarginals<-function(res.list, ## rest of arguments as for call to C fitabn
                       data.df,dag.m,var.types,max.parents,
                       mean,prec,loggam.shape,loggam.inv.scale,
                       max.iters,epsabs,verbose,error.verbose,
                       grouped.vars,## int.vector of variables which are mixed model nodes
                       group.ids,
                       epsabs.inner,max.iters.inner,finite.step.size,hessian.params,max.iters.hessian,min.pdf,marginal.node,marginal.param, variate.vec,n.grid){

if(!is.null(marginal.node) && is.null(variate.vec)){stop("must supply variate.vec if using a single node!");}

marginals<-list();

if( !is.null(marginal.node)) {## in single node case 
    if(res.list[["error.code"]][marginal.node]!=0){stop("---- Cannot compute marginal density as the mlik is unreliable for this node.\nTry re-running with different finite difference parameters ----");}
    nodeid<-marginal.node;
    paramid<-marginal.param;
    cat("processing ",names(res.list$modes[[nodeid]])[paramid],"\n");
    curnom<-names(res.list$modes[[nodeid]])[paramid];
    first<-TRUE;
    for(betafixed in variate.vec){
               marg.res <- .Call("fitabn_marginals",data.df,as.integer(dag.m),as.integer(dim(data.df)[2]),
                            as.integer(var.types),as.integer(max.parents),
                            as.double(mean),as.double(1/sqrt(prec)),as.double(loggam.shape),as.double(1/loggam.inv.scale),
                            as.integer(max.iters),as.double(epsabs),
                            as.integer(verbose),as.integer(error.verbose),
                            grouped.vars,## int.vector of variables which are mixed model nodes
                            group.ids,
                            as.double(epsabs.inner),
                            as.integer(max.iters.inner),
                            as.double(finite.step.size),
                            as.double(hessian.params),
                            ## additional marginal arguments - node, parameter, and relevant modes
                            as.integer(nodeid-1),## C index from 0
                            as.integer(paramid-1),## C index from 0
                            as.double(res.list$modes[[nodeid]]), ## relevant mode estimates
                            as.double(betafixed), ## value to be evaluated e.g. want f(x) given x - betafixed is x
                            as.double(res.list[[names(res.list$modes)[nodeid]]]),## mlik for the node
                            as.integer(max.iters.hessian),
              PACKAGE="abn" ## uncomment to load as package not shlib
              );
          if(first){tmp<-matrix(data=c(betafixed,marg.res),ncol=2,byrow=TRUE);colnames(tmp)<-c("x","f(x)");
                    marginals[[curnom]]<-tmp;first<-FALSE;
          } else{ marginals[[curnom]]<-rbind(marginals[[curnom]],c(betafixed,marg.res));}

          }
         marginals[[curnom]]<-marginals[[curnom]][order(marginals[[curnom]][,1]),];


} else {

  ## for each node
  for(nodeid in 1:length(res.list$modes)){
   #  nodeid<-2;
     ## for each parameter
     for(paramid in 1:length(res.list$modes[[nodeid]])){
   #  paramid<-2;    
         cat("processing ",names(res.list$modes[[nodeid]])[paramid],"\n");      
         curnom<-names(res.list$modes[[nodeid]])[paramid]; 
         if(res.list[["error.code"]][nodeid]!=0){cat("--- NOTE DROPPING parameter: ",curnom," as mlik is unreliable for this node.\nTry re-running with different finite difference parameters---\n");}  
         betafixedMode<-res.list$modes[[nodeid]][paramid];## just evaluate at the mode for paramid in node id
         ## STEP 1. get the value at the mode
         ## marg.res is a single number
         betafixed<-betafixedMode;
         marg.res <- .Call("fitabn_marginals",data.df,as.integer(dag.m),as.integer(dim(data.df)[2]),
                            as.integer(var.types),as.integer(max.parents),
                            as.double(mean),as.double(1/sqrt(prec)),as.double(loggam.shape),as.double(1/loggam.inv.scale),
                            as.integer(max.iters),as.double(epsabs),
                            as.integer(verbose),as.integer(error.verbose),
                            grouped.vars,## int.vector of variables which are mixed model nodes
                            group.ids,
                            as.double(epsabs.inner),
                            as.integer(max.iters.inner),
                            as.double(finite.step.size),
                            as.double(hessian.params),
                            ## additional marginal arguments - node, parameter, and relevant modes
                            as.integer(nodeid-1),## C index from 0
                            as.integer(paramid-1),## C index from 0
                            as.double(res.list$modes[[nodeid]]), ## relevant mode estimates
                            as.double(betafixed), ## value to be evaluated e.g. want f(x) given x - betafixed is x
                            as.double(res.list[[names(res.list$modes)[nodeid]]]),## mlik for the node
                            as.integer(max.iters.hessian),
              PACKAGE="abn" ## uncomment to load as package not shlib
              );

         gvalue.mode<-marg.res;
         tmp<-matrix(data=c(betafixedMode,gvalue.mode),ncol=2,byrow=TRUE);colnames(tmp)<-c("x","f(x)");
         marginals[[curnom]]<-tmp;
         #print(tmp);

         ## STEP 2. - try a grid spreading from the mode 
         if(betafixedMode<1E-05){## got a "zero mode 
                                try.these<-c(seq(betafixedMode,by= -0.05,length=3),seq(betafixedMode,by=0.05,length=3));
                                try.these<-try.these[order(try.these)];
         } else {try.these<-c(seq(betafixedMode,by= -betafixedMode*0.1,length=3),seq(betafixedMode,by=betafixedMode*0.1,length=3));
                                try.these<-try.these[order(try.these)];}

         for(betafixed in try.these){
               marg.res <- .Call("fitabn_marginals",data.df,as.integer(dag.m),as.integer(dim(data.df)[2]),
                            as.integer(var.types),as.integer(max.parents),
                            as.double(mean),as.double(1/sqrt(prec)),as.double(loggam.shape),as.double(1/loggam.inv.scale),
                            as.integer(max.iters),as.double(epsabs),
                            as.integer(verbose),as.integer(error.verbose),
                            grouped.vars,## int.vector of variables which are mixed model nodes
                            group.ids,
                            as.double(epsabs.inner),
                            as.integer(max.iters.inner),
                            as.double(finite.step.size),
                            as.double(hessian.params),
                            ## additional marginal arguments - node, parameter, and relevant modes
                            as.integer(nodeid-1),## C index from 0
                            as.integer(paramid-1),## C index from 0
                            as.double(res.list$modes[[nodeid]]), ## relevant mode estimates
                            as.double(betafixed), ## value to be evaluated e.g. want f(x) given x - betafixed is x
                            as.double(res.list[[names(res.list$modes)[nodeid]]]),## mlik for the node
                            as.integer(max.iters.hessian),
              PACKAGE="abn" ## uncomment to load as package not shlib
              );
          marginals[[curnom]]<-rbind(marginals[[curnom]],c(betafixed,marg.res));

          }
         marginals[[curnom]]<-marginals[[curnom]][order(marginals[[curnom]][,1]),];

     ### have we found a suitably wide range
     gleft<-marginals[[curnom]][1,2];
     gright<-marginals[[curnom]][dim(marginals[[curnom]])[1],2];
    #cat("gleft=",gleft,"gright=",gright,"\n");
  
   ### STEP 3. it may be that the initial grid is to crude resulting in just a "spike" so make this grid finer according to the following rule
   max.val.index<-which(marginals[[curnom]][,2]==max(marginals[[curnom]][,2]));

   if(min(abs(marginals[[curnom]][-max.val.index,2]-gvalue.mode))>0.5*gvalue.mode ){## this is the rule 

    if(betafixedMode<1E-05){## got a "zero mode 
                                try.these<-c(seq(betafixedMode,by= -0.004,length=3),seq(betafixedMode,by=0.004,length=3));
                                try.these<-try.these[order(try.these)];
         } else {try.these<-c(seq(betafixedMode,by= -betafixedMode*0.004,length=3),seq(betafixedMode,by=betafixedMode*0.004,length=3));
                                try.these<-try.these[order(try.these)];}

         for(betafixed in try.these){
               marg.res <- .Call("fitabn_marginals",data.df,as.integer(dag.m),as.integer(dim(data.df)[2]),
                            as.integer(var.types),as.integer(max.parents),
                            as.double(mean),as.double(1/sqrt(prec)),as.double(loggam.shape),as.double(1/loggam.inv.scale),
                            as.integer(max.iters),as.double(epsabs),
                            as.integer(verbose),as.integer(error.verbose),
                            grouped.vars,## int.vector of variables which are mixed model nodes
                            group.ids,
                            as.double(epsabs.inner),
                            as.integer(max.iters.inner),
                            as.double(finite.step.size),
                            as.double(hessian.params),
                            ## additional marginal arguments - node, parameter, and relevant modes
                            as.integer(nodeid-1),## C index from 0
                            as.integer(paramid-1),## C index from 0
                            as.double(res.list$modes[[nodeid]]), ## relevant mode estimates
                            as.double(betafixed), ## value to be evaluated e.g. want f(x) given x - betafixed is x
                            as.double(res.list[[names(res.list$modes)[nodeid]]]),## mlik for the node
                            as.integer(max.iters.hessian),
              PACKAGE="abn" ## uncomment to load as package not shlib
              );
          marginals[[curnom]]<-rbind(marginals[[curnom]],c(betafixed,marg.res));

          }
         marginals[[curnom]]<-marginals[[curnom]][order(marginals[[curnom]][,1]),];

        } ## end of finer grid

       ### STEP 4. compute an increment based on a simple linear slope estimate (separate for left and right directions and which is roughly
       ###         equivalent to a drop of 10% in height for each step across the grid
             left<-try.these[1];
             right<-try.these[length(try.these)];
             slope.to.left<-(gvalue.mode-gleft)/abs(betafixedMode-left);
             slope.to.right<-(gvalue.mode-gright)/abs(betafixedMode-right);
             
             ## what stepsize is equivalent to a 10% drop in height
             delta.left<-(0.1*gvalue.mode)/slope.to.left;
             delta.right<-(0.1*gvalue.mode)/slope.to.right;
             gvalue<-.Machine$double.xmax;
             betafixed<-betafixedMode;
     #cat("delta.left=",delta.left," delta.right=",delta.right,"\n");
      ### STEP 5. perform the iteration to the left       
             while(gvalue>min.pdf){

             betafixed<- betafixed - delta.left;
             marg.res <- .Call("fitabn_marginals",data.df,as.integer(dag.m),as.integer(dim(data.df)[2]),
                            as.integer(var.types),as.integer(max.parents),
                            as.double(mean),as.double(1/sqrt(prec)),as.double(loggam.shape),as.double(1/loggam.inv.scale),
                            as.integer(max.iters),as.double(epsabs),
                            as.integer(verbose),as.integer(error.verbose),
                            grouped.vars,## int.vector of variables which are mixed model nodes
                            group.ids,
                            as.double(epsabs.inner),
                            as.integer(max.iters.inner),
                            as.double(finite.step.size),
                            as.double(hessian.params),
                            ## additional marginal arguments - node, parameter, and relevant modes
                            as.integer(nodeid-1),## C index from 0
                            as.integer(paramid-1),## C index from 0
                            as.double(res.list$modes[[nodeid]]), ## relevant mode estimates
                            as.double(betafixed), ## value to be evaluated e.g. want f(x) given x - betafixed is x
                            as.double(res.list[[names(res.list$modes)[nodeid]]]),## mlik for the node
                            as.integer(max.iters.hessian),
              PACKAGE="abn" ## uncomment to load as package not shlib
              );
          gvalue<-marg.res;
          marginals[[curnom]]<-rbind(marginals[[curnom]],c(betafixed,marg.res));
          }

       ### STEP 6. perform the iteration to the right   
            gvalue<-.Machine$double.xmax;
             betafixed<-betafixedMode;
             
             while(gvalue>min.pdf){

             betafixed<- betafixed + delta.right;
             marg.res <- .Call("fitabn_marginals",data.df,as.integer(dag.m),as.integer(dim(data.df)[2]),
                            as.integer(var.types),as.integer(max.parents),
                            as.double(mean),as.double(1/sqrt(prec)),as.double(loggam.shape),as.double(1/loggam.inv.scale),
                            as.integer(max.iters),as.double(epsabs),
                            as.integer(verbose),as.integer(error.verbose),
                            grouped.vars,## int.vector of variables which are mixed model nodes
                            group.ids,
                            as.double(epsabs.inner),
                            as.integer(max.iters.inner),
                            as.double(finite.step.size),
                            as.double(hessian.params),
                            ## additional marginal arguments - node, parameter, and relevant modes
                            as.integer(nodeid-1),## C index from 0
                            as.integer(paramid-1),## C index from 0
                            as.double(res.list$modes[[nodeid]]), ## relevant mode estimates
                            as.double(betafixed), ## value to be evaluated e.g. want f(x) given x - betafixed is x
                            as.double(res.list[[names(res.list$modes)[nodeid]]]),## mlik for the node
                            as.integer(max.iters.hessian),
              PACKAGE="abn" ## uncomment to load as package not shlib
              );
          gvalue<-marg.res;
          marginals[[curnom]]<-rbind(marginals[[curnom]],c(betafixed,marg.res));
          }

          marginals[[curnom]]<-marginals[[curnom]][order(marginals[[curnom]][,1]),];   
          
    ## STEP 7. now drop the tails in which the g value is below the min
          keep.these<-which(marginals[[curnom]][,2]>=min.pdf);
          keep.these<-c(min(keep.these)-1,keep.these,max(keep.these)+1);## adjust so take next variate beyond cut-off
          marginals[[curnom]]<-marginals[[curnom]][keep.these,];
 
    ## STEP 8. if asked for go this on an equal grid of n.grid points

   
   if(!is.null(n.grid)){
   mymarg<-matrix(data=rep(NA,n.grid*2),ncol=2,byrow=TRUE);colnames(mymarg)<-c("x","f(x)");
     ## repeat entirely using equally spaced grid from min x to max x
        i<-1; 
        for(betafixed in seq(from=marginals[[curnom]][1,1],to=marginals[[curnom]][dim(marginals[[curnom]])[1],1],len=n.grid)){
               marg.res <- .Call("fitabn_marginals",data.df,as.integer(dag.m),as.integer(dim(data.df)[2]),
                            as.integer(var.types),as.integer(max.parents),
                            as.double(mean),as.double(1/sqrt(prec)),as.double(loggam.shape),as.double(1/loggam.inv.scale),
                            as.integer(max.iters),as.double(epsabs),
                            as.integer(verbose),as.integer(error.verbose),
                            grouped.vars,## int.vector of variables which are mixed model nodes
                            group.ids,
                            as.double(epsabs.inner),
                            as.integer(max.iters.inner),
                            as.double(finite.step.size),
                            as.double(hessian.params),
                            ## additional marginal arguments - node, parameter, and relevant modes
                            as.integer(nodeid-1),## C index from 0
                            as.integer(paramid-1),## C index from 0
                            as.double(res.list$modes[[nodeid]]), ## relevant mode estimates
                            as.double(betafixed), ## value to be evaluated e.g. want f(x) given x - betafixed is x
                            as.double(res.list[[names(res.list$modes)[nodeid]]]),## mlik for the node
                            as.integer(max.iters.hessian),
              PACKAGE="abn" ## uncomment to load as package not shlib
              );
           mymarg[i,]<-c(betafixed,marg.res);
           i<-i+1;
             }
          ## overwrite existing marginal computation....
          marginals[[curnom]]<-mymarg;
                }

             } #parameter

        } #node

     } ## end of if single node or all

     #cat("f(",betafixed,")=",marg.res,"\n");   
     

res.list[["marginals"]]<-marginals;


return(res.list);


}
                                                
