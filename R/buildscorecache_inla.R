###############################################################################
## fitabn.R --- 
## Author          : Fraser Lewis
## Last modified   : 03/08/2012
###############################################################################

## fit a given DAG to data using the distribution given
buildscorecache.inla <- function(data.df=NULL, data.dists=NULL,ntrials=NULL, exposure=NULL, group.var=NULL,cor.vars=NULL,
                            dag.banned=NULL, dag.retained=NULL,max.parents=NULL,which.nodes=NULL,
                            defn.res=NULL,dry.run=FALSE,
                            verbose=FALSE,centre=TRUE,mean=0, prec=0.001,loggam.shape=1,loggam.inv.scale=5e-05){
      
      use.inla<-TRUE;
      ## check grouping variables
      list.group.var<-check.groups(group.var,data.df,cor.vars,use.inla);## returns ammended data.df and suitable variables
      data.df<-list.group.var$data.df;## this has removes the grouping variable from data.df
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
      max.parents<-check.parents(data.df,max.parents,use.inla,group.var);## returns an integer vector of the same length of number of nodes

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
      ## use inla to compute scores
      #########################################################
         if(!require(INLA)){stop("library INLA is not available!\nR-INLA is available from http://www.r-inla.org/download\nAfter installation please use inla.upgrade() to get the latest version (this is required)");}
       
       dag.m.loc<-dag.banned;## create a dummy dag since calc.node.mlik expects a full DAG not just a row
       dag.m.loc[,]<-0;## fills dag.m with zeros
       ## set up other variables needed for inla() but which are fixed here
       compute.fixed.loc<-FALSE;
       mean.intercept.loc<-mean; prec.intercept.loc<-prec; mean.loc<-mean; prec.loc<-prec; 
       loggam.shape.loc<-loggam.shape;
       loggam.inv.scale.loc<-loggam.inv.scale;
       scores<-rep(NA,length(defn.res[[1]]));
       
        node.num<-1;
        for(child in defn.res[["children"]]){ ## each entry is a node ID (1 through numNodes)
                   cat("processing...",node.num," of ",length(defn.res[[1]]),"\n",sep="");
                   #child.name<-colnames(defn.res[["node.defn"]])[child];        
                   dag.m.loc[,]<-0;## reset to empty DAG
                   dag.m.loc[child,]<-defn.res[["node.defn"]][node.num,];## copy node define into DAG  

                   if(child%in%(grouped.vars+1)){#cat("got grouped node - need mixed model\n");## got a mixed node
                            scores[node.num]<-calc.node.mlik.inla.mixed(group.var,
                                                                        child.loc=child,
                                                                        dag.m.loc=dag.m.loc,
                                                                        data.df.loc=data.df,
                                                                        data.dists.loc=data.dists,
                                                                        ntrials.loc=mylist$ntrials,
                                                                        exposure.loc=mylist$exposure, 
                                                                        compute.fixed.loc=compute.fixed.loc,
                                                                        mean.intercept.loc= mean.intercept.loc, 
                                                                        prec.intercept.loc=prec.intercept.loc, 
                                                                        mean.loc=mean.loc,
                                                                        prec.loc=prec.loc,
                                                                        loggam.shape.loc=loggam.shape.loc,
                                                                        loggam.inv.scale.loc=loggam.inv.scale.loc,
                                                                        verbose.loc=verbose);

                   } else {scores[node.num]<-calc.node.mlik.inla(child.loc=child,
                                                          dag.m.loc=dag.m.loc,
                                                          data.df.loc=data.df,
                                                          data.dists.loc=data.dists,
                                                          ntrials.loc=mylist$ntrials,
                                                          exposure.loc=mylist$exposure, 
                                                          compute.fixed.loc=compute.fixed.loc,
                                                          mean.intercept.loc=mean.intercept.loc,
                                                          prec.intercept.loc=prec.intercept.loc,
                                                          mean.loc=mean.loc,
                                                          prec.loc=prec.loc,
                                                          loggam.shape.loc=loggam.shape.loc,
                                                          loggam.inv.scale.loc=loggam.inv.scale.loc,
                                                          verbose.loc=verbose);
               }
                   node.num<-node.num+1;## update index in complete node list
      }

                  defn.res[["mlik"]]<-scores;
                  #defn.res[["banned.arcs"]]<-dag.banned;
                  #defn.res[["retained.arcs"]]<-dag.retained;
                
                  ## final part as we want a version of data.df without the grouping variable for use with mostprobable() etc
                  use.inla<-FALSE;
                  list.group.var<-check.groups(group.var,data.df,cor.vars,use.inla);## returns ammended data.df and suitable variables
                  data.df<-list.group.var$data.df;## this has removes the grouping variable from data.df
                  defn.res[["data.df"]]<-data.df;
                 
     


return(defn.res);

}
