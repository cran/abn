###############################################################################
## fitabn.R --- 
## Author          : Fraser Lewis
## Last modified   : 29/08/2012
###############################################################################

## fit a given DAG to data using the distribution given
fitabn.inla <- function(dag.m=NULL, data.df=NULL, data.dists=NULL, group.var=NULL,cor.vars=NULL,create.graph=FALSE,compute.fixed=FALSE,
                   mean=0, prec=0.001,loggam.shape=1,loggam.inv.scale=5e-05,verbose=FALSE, centre=TRUE)
                   {
      ntrials<-NULL;
      exposure<-NULL;
      use.inla<-TRUE;
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
      ## use inla
      #########################################################
      
       if(!require(INLA)){stop("library INLA is not available!\nR-INLA is available from http://www.r-inla.org/download\nAfter installation please use inla.upgrade() to get the latest version (this is required)");}
       mean.intercept<-mean;## use same as for rest of linear terms 
       prec.intercept<-prec;## use same as for rest of linear terms
       for(child in 1:dim(dag.m)[1]){ ## for each node in network fit the appropriate model using call in inla
                   child.name<-colnames(dag.m)[child]; 
                   if(child%in%(grouped.vars+1)){#cat("got grouped node - need mixed model\n");## got a mixed node
                            res.list[[child.name]]<-calc.node.mlik.inla.mixed(group.var,child,dag.m,data.df,data.dists,mylist$ntrials,mylist$exposure, compute.fixed, mean.intercept, 
                                                               prec.intercept, mean, prec,loggam.shape,loggam.inv.scale,verbose);
                   } else {## just need glm 
                   res.list[[child.name]]<-calc.node.mlik.inla(child,dag.m,data.df,data.dists,mylist$ntrials,mylist$exposure, compute.fixed, mean.intercept, 
                                                               prec.intercept, mean, prec,loggam.shape,loggam.inv.scale,verbose);}
       }

      ## Now get marginal distributions if requested. Complicated by the fact that re-organization of the inla output is required 

      if(!compute.fixed){res.list[["mlik"]]<-sum(unlist(res.list));## add in total mlik for DAG
                        # res.list[["use.inla"]]<-TRUE;
      } else {
             ## We want marginals - which have already been computed above so reorganise this into format consistent with C output (below)
             newlist<-list();
             for(i in colnames(dag.m)){newlist[[i]]<-res.list[[i]]$mlik[2];} ## [2] is for the gaussian mlik est. from INLA
             newlist[["mlik"]]<-sum(unlist(newlist));
             inner.list<-list();nom<-NULL
             for(i in colnames(dag.m)){## for each covariate
                 #names(res.list[[i]]$marginals.fixed)<-paste(i,"|",names(res.list[[i]]$marginals.fixed),sep="");
                 nom<-c(nom,paste(i,"|",names(res.list[[i]]$marginals.fixed),sep=""));
                 if(data.dists[[i]]=="gaussian" && !(i%in%cor.vars)){## got gaussian variable and it is not grouped
                                        #names(res.list[[i]]$marginals.hyperpar)<-paste(i,"|precision",sep="");
                                         nom<-c(nom,paste(i,"|precision",sep=""));
                                         }
                 if(data.dists[[i]]=="gaussian" && i%in%cor.vars){## the grouped variable 
                                         nom<-c(nom,paste(i,"|precision",sep=""));
                                         nom<-c(nom,paste(i,"|Groups-precision",sep=""));
                            }

                 if(data.dists[[i]]=="binomial" && i%in%cor.vars){## the grouped variable 
                                         nom<-c(nom,paste(i,"|Groups-precision",sep=""));
                            } 
                 if(data.dists[[i]]=="poisson" && i%in%cor.vars){## the grouped variable 
                                         nom<-c(nom,paste(i,"|Groups-precision",sep=""));
                            }
              }
             index=1;
             for(i in colnames(dag.m)){## for each variable
                 tmp.res<-res.list[[i]]$marginals.fixed;
                    for(j in 1:length(tmp.res)){inner.list[index]<-tmp.res[j]; index<-index+1;
                                                }
                    if(data.dists[[i]]=="gaussian" && !(i%in%cor.vars)){## got gaussian variable and it is not grouped
                                        inner.list[index]<-res.list[[i]]$marginals.hyperpar[1];index<-index+1;}

                    if(data.dists[[i]]=="gaussian" && i%in%cor.vars){## the grouped variable 
                                        inner.list[index]<-res.list[[i]]$marginals.hyperpar[1];index<-index+1;## the residual precision
                                        inner.list[index]<-res.list[[i]]$marginals.hyperpar[2];index<-index+1;## group level precision
                                         }
       
                    if(data.dists[[i]]=="binomial" && i%in%cor.vars){## the grouped variable 
                                        inner.list[index]<-res.list[[i]]$marginals.hyperpar[1];index<-index+1;}
                    
                    if(data.dists[[i]]=="poisson" && i%in%cor.vars){## the grouped variable 
                                        inner.list[index]<-res.list[[i]]$marginals.hyperpar[1];index<-index+1;} 
                   
             }
             #newlist[["use.inla"]]<-TRUE;
             names(inner.list)<-nom;
             newlist[["marginals"]]<-inner.list;
             for(i in 1:length(newlist$marginals)){colnames(newlist$marginals[[i]])<-c("x","pdf(x)");}
             #return(inner.list);
             res.list<-newlist;## overwrite
             #rm(res.list);
             #return(newlist);
              }

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
