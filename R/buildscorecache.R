###############################################################################
## fitabn.R --- 
## Author          : Fraser Lewis
## Last modified   : 03/08/2012
###############################################################################

## fit a given DAG to data using the distribution given
buildscorecache <- function(data.df=NULL, data.dists=NULL, group.var=NULL,cor.vars=NULL,
                            dag.banned=NULL, dag.retained=NULL,max.parents=NULL,
                            which.nodes=NULL,defn.res=NULL,dry.run=FALSE,
                            verbose=FALSE,centre=TRUE,mean=0, prec=0.001,loggam.shape=1,
                            loggam.inv.scale=5e-05, max.iters=100,
                            epsabs=1e-7,error.verbose=FALSE,term.output.freq=100,
                            epsabs.inner=1e-6,max.iters.inner=100,
                            finite.step.size=1e-7,
                            hessian.params=c(1E-04,1E-02),max.iters.hessian=10){
      
      use.inla<-FALSE;
      ntrials<-NULL;
      exposure<-NULL;
      ## check grouping variables
      list.group.var<-check.groups(group.var,data.df,cor.vars,use.inla);## returns ammended data.df and suitable variables
      data.df<-list.group.var$data.df;## this has removed the grouping variable from data.df
      grouped.vars<-list.group.var$grouped.vars;## int vect of variables to be treated as grouped
      group.ids<-list.group.var$group.ids;## int vector of group membership ids
      #data.df.inla<-list.group.var$data.df.inla;

      ### run a series of checks on the data and distributions passed
      mylist<-check.data(data.df,data.dists,ntrials,exposure,use.inla,group.var);## return a list with entries bin, gaus, pois, ntrials and exposure

      ## run a series of common sense checks on the banned DAG and retain DAG
      if(!is.null(dag.banned)){check.dag(dag.banned,data.df=data.df,is.ban.matrix=TRUE,use.inla,group.var);
      } else {dag.banned<-check.dag(dag.banned,data.df=data.df,is.ban.matrix=TRUE,use.inla,group.var);} ##if null just create empty mat and return

      if(!is.null(dag.retained)){check.dag(dag.retained,data.df=data.df,is.ban.matrix=FALSE,use.inla,group.var);
      } else {dag.retained<-check.dag(dag.retained,data.df=data.df,is.ban.matrix=FALSE,use.inla,group.var);} ##if null just create empty mat and return


      ## check max.parents is a list with suitable entries      
      if(is.null(defn.res)){## not supplying custom parent sets
        max.parents<-check.parents(data.df,max.parents,use.inla,group.var);## returns an integer vector of the same length of number of nodes
      } else {## providing custom parent set and so max.parents is irrelevant
              max.parents<-max(apply(defn.res[["node.defn"]],1,sum));}## note this is just a dummy and not actually used anywhere!
      ## check retain does not ask for more arcs to be retained than allowed in max.parents
      max.retain<-apply(dag.retained,1,sum);## number of parents per node to retain
      if(length(which( (max.retain>max.parents) == TRUE))>0){stop("dag.retained is inconsistent with max.parents!");}

      ## check that arcs than are banned are also not retained      
      if(length(which(which(as.integer(dag.banned)==1)%in%which(as.integer(dag.retained)==1)==TRUE))>0){stop("dag.banned and dag.retained are inconsistent!");}
      
      ## check which.nodes is sensible
      if(is.null(defn.res)){which.nodes<-check.which.nodes(data.df,which.nodes,use.inla,group.var);
      } else {## have user supplied children and parent combinations 
              which.nodes<-unique(defn.res$children);}

      ## after checks do some coercion and data changes
      ## coerce binary factors to become 0/1 integers - the 0 is based on the first entry in levels()
      if(!is.null(mylist$bin)){## have at least one binary variable
        for(i in mylist$bin){data.df[,i]<-as.numeric(data.df[,i])-1;}
      }

      ## standardize gaussian variables to zero mean and sd=1
      if(centre && !is.null(mylist$gaus)){## have at least one binary variable
        for(i in mylist$gaus){data.df[,i]<-(data.df[,i]-mean(data.df[,i]))/sd(data.df[,i]);}
      }

      ## down to here we have all the data correct and now call C buildnodecache() to create all the node definitions. There is a separate R buildnodecache()
      ## function which we have just replicated here - rather than call it explicitly since we then has to faff about passing arguments etc.
      if(is.null(defn.res)){
      ## pass to C the number (number_of_nodes,banned_arc_as_vector, retain_arcs_as_vector, max_parents_as_vector
      res<-.Call("buildcachematrix",dim(dag.banned)[1],as.integer(dag.banned),as.integer(dag.retained), max.parents, which.nodes
              ,PACKAGE="abn" ## uncomment to load as package not shlib
              ) 
      defn.res<-list();
      defn.res[["children"]]<-res[[1]];
      if(use.inla && !is.null(group.var) ){
      defn.res[["node.defn"]]<-matrix(data=res[[2]],byrow=TRUE,ncol=dim(data.df)[2]-1);## because data.df has one too many columns
      colnames(defn.res[["node.defn"]])<-names(data.df)[-dim(data.df)[2]];## drop last col since this is the grouping var
      } else { defn.res[["node.defn"]]<-matrix(data=res[[2]],byrow=TRUE,ncol=dim(data.df)[2]);
               colnames(defn.res[["node.defn"]])<-names(data.df);}
      rm(res);
      
      } else { ## some check since user has supplied defn.res
               if(!is.list(defn.res)){stop("defn.res must be a list");}
               if(length(defn.res)!=2){stop("defn.res must have two entries");}
               if(!(is.vector(defn.res[[1]]) && is.matrix(defn.res[[2]]))){stop("defn.res is wrong format");}
               if(!(max(defn.res[[2]])==1 && min(defn.res[[2]])==0)){stop("defn.res is wrong format - must only be 0,1 in node definitions");} 
              }
      
      if(dry.run){## don't do any computation just return the node definitions
                  cat("No computation - returning only the node combinations\n"); 
                  return(defn.res);}
      #compute.fixed=FALSE,
      #             mean.intercept=0, prec.intercept=0.001, mean=0, prec=0.001,verbose=FALSE, centre=TRUE
      
      ################################################
      ## we now have a set of models - single node - definitions and want each one evaluated using a call to inla() or internal C code
      ## loop through each node ...etc
      #################################################

      #########################################################
      ## use C to compute scores
      #########################################################
      cache.defn<-as.integer(defn.res[["node.defn"]]);## into one long vector filled by col
      children<-as.integer(defn.res[["children"]]);
      numVars<-as.integer(dim(defn.res[["node.defn"]])[2]);## number of variables in the DAG
      numRows<-as.integer(dim(defn.res[["node.defn"]])[1]);## total number of different node-parent combinations
      numparents.per.node<-as.integer(table(children));## number of parent combinations per variable
      #return(numparents.per.node);

      for(i in 1:dim(data.df)[2]){data.df[,i]<-as.double(data.df[,i]);}#coerce ALL cols to double equivalents
      var.types<-get.var.types(data.dists); ## get distributions in terms of a code
      max.parents<-max(apply(defn.res[[2]],1,sum));
      res <- .Call("fitnodes",data.df,children,cache.defn,numVars,numRows,numparents.per.node,
                            as.integer(var.types),as.integer(max.parents),
                            as.double(mean),as.double(1/sqrt(prec)),as.double(loggam.shape),as.double(1/loggam.inv.scale),
                                         as.integer(max.iters),as.double(epsabs),
                                         as.integer(verbose),as.integer(error.verbose),
                             which.nodes, as.integer(term.output.freq),
                             grouped.vars,## new stuff from here 30/09/2012. int.vector of variables which are mixed model nodes
                             group.ids,
                             as.double(epsabs.inner),
                             as.integer(max.iters.inner),
                             as.double(finite.step.size),
                             as.double(hessian.params),
                             as.integer(max.iters.hessian)

                   ,PACKAGE="abn" ## uncomment to load as package not shlib
                  )
      #return(res);
      defn.res[["mlik"]]<-res[[1]];
      defn.res[["error.code"]]<-res[[2]];
      defn.res[["error.code.desc"]]<-as.character(defn.res[["error.code"]]);
      defn.res[["error.code.desc"]]<-ifelse(defn.res[["error.code.desc"]]==0,"success",defn.res[["error.code.desc"]]);
      defn.res[["error.code.desc"]]<-ifelse(defn.res[["error.code.desc"]]==1,"warning: mode results may be unreliable (optimiser terminated unusually)",defn.res[["error.code.desc"]]);
      defn.res[["error.code.desc"]]<-ifelse(defn.res[["error.code.desc"]]==2,"warning: initial interval in fin.diff unusable - h_guess used",defn.res[["error.code.desc"]]);
      defn.res[["error.code.desc"]]<-ifelse(defn.res[["error.code.desc"]]==3,"warning: fin.diff adaptive step size failed (hessian NaN or Inf) - h_guess used",defn.res[["error.code.desc"]]);
      defn.res[["error.code.desc"]]<-ifelse(defn.res[["error.code.desc"]]==4,"warning: fin.diff terminated unusually - h_guess used",defn.res[["error.code.desc"]]);

      defn.res[["hessian.accuracy"]]<-res[[3]];
      #defn.res[["banned.arcs"]]<-dag.banned;
      #defn.res[["retained.arcs"]]<-dag.retained;

      defn.res[["data.df"]]<-data.df;
      


return(defn.res);

}
