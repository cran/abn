### R code from vignette source 'abn_v0.85.Rnw'

###################################################
### code chunk number 1: abn_v0.85.Rnw:113-151
###################################################
library( abn) # load library
bin.nodes<-c( 1,3,4,6,9,10,11,12,15,18,19,20,21,26,27,28,32); 
var33.cat<-var33[,bin.nodes]; #categorical nodes only

mydag<-matrix(c(  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v1
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v3
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v4  
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v6  
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v9  
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v10  
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v11  
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v12  
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v15  
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v18 
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v19
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v20 
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v21 
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v26 
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v27 
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #v28 
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0  #v32 
               ),byrow=TRUE,ncol=17); 
colnames( mydag)<-rownames( mydag)<-names( var33.cat);#set names
## move back to independence model
mydag["v11","v12"]<-0;mydag["v11","v10"]<-0;mydag["v4","v3"]<-0;
## setup distribution list for each categorical node
mydists.cat<-list( v1 ="binomial", v3 = "binomial",
       v4 = "binomial",  v6 = "binomial",  v9 = "binomial",
      v10 = "binomial", v11 = "binomial", v12 = "binomial",
      v15 = "binomial", v18 = "binomial", v19 = "binomial",
      v20 = "binomial", v21 = "binomial", v26 = "binomial",
      v27 = "binomial", v28 = "binomial", v32 = "binomial");
ind.mod.cat <- fitabn( data.df=var33.cat, dag.m=mydag,
                      data.dists=mydists.cat, verbose=FALSE); 
## change to verbose=TRUE if one want to check how change the 
## score for each individual node
ind.mod.cat$mlik 
## network score for a model with conditional independencies


###################################################
### code chunk number 2: abn_v0.85.Rnw:154-160
###################################################
## now fit the model with some conditional dependencies 
mydag["v11","v12"]<-1;mydag["v11","v10"]<-1;mydag["v4","v3"]<-1;
dep.mod.cat <- fitabn( data.df=var33.cat, dag.m=mydag,
                      data.dists=mydists.cat, verbose=FALSE);
dep.mod.cat$mlik
## network score for a model with conditional dependencies


###################################################
### code chunk number 3: abn_v0.85.Rnw:163-168
###################################################
tographviz( dag=mydag, data.df=var33.cat, data.dists=mydists.cat,
           outfile="mydagcat.dot", directed=TRUE);#create file
# mydagcat.dot can then be processed with graphviz
# unix shell "dot -Tpdf mydagcat.dot -o mydagcat.pdf" 
# or use gedit if on Windows


###################################################
### code chunk number 4: abn_v0.85.Rnw:181-197
###################################################
var33.cts<-var33[,-bin.nodes]; #drop categorical nodes
mydag<-matrix( 0, 16, 16);
colnames( mydag)<-rownames( mydag)<-names( var33.cts);#set names
## setup distribution list for each continuous node
mydists.cts<-list( v2 = "gaussian", v5 = "gaussian",
        v7 = "gaussian", v8 = "gaussian", v13 = "gaussian", 
        v14 = "gaussian", v16 = "gaussian", v17 = "gaussian",
        v22 = "gaussian", v23 = "gaussian", v24 = "gaussian",
        v25 = "gaussian", v29 = "gaussian", v30 = "gaussian",
        v31 = "gaussian", v33 = "gaussian");
## now fit the model defined in mydag - full independence
ind.mod.cts <- fitabn( data.df=var33.cts, dag.m=mydag,
               data.dists=mydists.cts, verbose=FALSE);
## uses default priors: N(mu=0,var=1000), 1/var=Gamma(0.001,1/0.001)
ind.mod.cts$mlik
# this is the network score=goodness of fit=log marginal likelihood


###################################################
### code chunk number 5: abn_v0.85.Rnw:200-214
###################################################
# now fit model with some conditional dependencies let v33 
## depend on v31, and v24 depend on 23, and v14 depend on v13
mydag["v33","v31"]<-1;
mydag["v24","v23"]<-1;
mydag["v14","v13"]<-1;
dep.mod.cts <- fitabn( data.df=var33.cts, dag.m=mydag,
               data.dists=mydists.cts, verbose=FALSE);
dep.mod.cts$mlik
# network score for a model with conditional independence
tographviz( dag=mydag, data.df=var33.cts, data.dists=mydists.cts,
            outfile="mydagcts.dot", directed=TRUE);#create file
# mydag.dot can then be processed with graphviz
# unix shell "dot -Tpdf mydagcts.dot -o mydagcts.pdf" or 
# use gedit if on Windows


###################################################
### code chunk number 6: abn_v0.85.Rnw:225-245
###################################################
mydag<-matrix( 0, 33, 33); 
colnames( mydag)<-rownames( mydag)<-names( var33);#set names
## setup distribution list for each mixed node
mydists.mix<-list( v1 = "binomial", v2 = "gaussian",
        v3 = "binomial", v4 = "binomial", v5 = "gaussian",
        v6 = "binomial", v7 = "gaussian", v8 = "gaussian",
        v9 = "binomial", v10 = "binomial", v11 = "binomial",
        v12 = "binomial", v13 = "gaussian", v14 = "gaussian",
        v15 = "binomial", v16 = "gaussian", v17 = "gaussian",
        v18 = "binomial", v19 = "binomial", v20 = "binomial",
        v21 = "binomial", v22 = "gaussian", v23 = "gaussian",
        v24 = "gaussian", v25 = "gaussian", v26 = "binomial",
        v27 = "binomial", v28 = "binomial", v29 = "gaussian",
        v30 = "gaussian", v31 = "gaussian", v32 = "binomial",
        v33 = "gaussian");
## now fit the model defined in mydag - full independence
ind.mod <- fitabn( data.df=var33, dag.m=mydag,
          data.dists=mydists.mix, verbose=FALSE);
ind.mod$mlik
# this is the network score with no conditional dependencies


###################################################
### code chunk number 7: abn_v0.85.Rnw:248-282
###################################################
# define a model with many independencies
mydag[2,1]<-1;
mydag[4,3]<-1;
mydag[6,4]<-1; mydag[6,7]<-1;
mydag[5,6]<-1;
mydag[7,8]<-1;  
mydag[8,9]<-1;
mydag[9,10]<-1;
mydag[11,10]<-1; mydag[11,12]<-1; mydag[11,19]<-1;
mydag[14,13]<-1;
mydag[17,16]<-1;mydag[17,20]<-1;
mydag[15,14]<-1; mydag[15,21]<-1;
mydag[18,20]<-1;
mydag[19,20]<-1;
mydag[21,20]<-1;
mydag[22,21]<-1;
mydag[23,21]<-1;
mydag[24,23]<-1;
mydag[25,23]<-1; mydag[25,26]<-1;
mydag[26,20]<-1;
mydag[33,31]<-1;
mydag[33,31]<-1;
mydag[32,21]<-1; mydag[32,31]<-1;mydag[32,29]<-1;    
mydag[30,29]<-1;
mydag[28,27]<-1; mydag[28,29]<-1;mydag[28,31]<-1;       
dep.mod <- fitabn( data.df=var33, dag.m=mydag, 
           data.dists=mydists.mix, verbose=FALSE);
dep.mod$mlik
# network score for a model with conditional independence
tographviz( dag=mydag, data.df=var33, data.dists=mydists.mix,
            outfile="mydag_all.dot", directed=TRUE);#create file
# mydag.dot can then be processed with graphviz
# unix shell "dot -Tpdf mydag_all.dot -o mydag_all.pdf" or use 
# gedit if on Windows


###################################################
### code chunk number 8: abn_v0.85.Rnw:302-323
###################################################
bin.nodes<-c(1,3,4,6,9,10,11,12,15,18,19,20,21,26,27,28,32); 
var33.cat<-var33[,bin.nodes];#categorical nodes only
mydag<-matrix( 0, 17, 17); 
colnames(mydag)<-rownames(mydag)<-names(var33.cat);#set names
## create banned and retain empty DAGs
banned.cat<-matrix( 0, 17, 17);
colnames(banned.cat)<-rownames(banned.cat)<-names(var33.cat);
retain.cat<-matrix( 0, 17, 17);
colnames(retain.cat)<-rownames(retain.cat)<-names(var33.cat);
## setup distribution list for each categorical node
mydists.cat<-list( v1 ="binomial", v3 = "binomial",
       v4 = "binomial",  v6 = "binomial",  v9 = "binomial",
      v10 = "binomial", v11 = "binomial", v12 = "binomial",
      v15 = "binomial", v18 = "binomial", v19 = "binomial",
      v20 = "binomial", v21 = "binomial", v26 = "binomial",
      v27 = "binomial", v28 = "binomial", v32 = "binomial");
## build cache of all the local computations
## this information is needed later when running a model search
mycache.cat<-buildscorecache( data.df=var33.cat,
             data.dists=mydists.cat, dag.banned=banned.cat, 
             dag.retained=retain.cat, max.parents=1);


###################################################
### code chunk number 9: abn_v0.85.Rnw:328-334 (eval = FALSE)
###################################################
## # Run a single search heuristic for an additive BN
## heur.res.cat<-search.hillclimber( score.cache=mycache.cat,
##               num.searches=1, seed=0, verbose=FALSE,
##               trace=FALSE, timing.on=FALSE);
## # Setting trace=TRUE, the majority consensus network is 
## # plotted as the searches progress


###################################################
### code chunk number 10: abn_v0.85.Rnw:339-358
###################################################
var33.cts<-var33[,-bin.nodes];#drop categorical nodes
mydag<-matrix( 0, 16, 16); 
colnames(mydag)<-rownames(mydag)<-names(var33.cts);#set names
banned.cts<-matrix( 0, 16, 16);
colnames(banned.cts)<-rownames(banned.cts)<-names(var33.cts);
retain.cts<-matrix( 0, 16, 16);
colnames(retain.cts)<-rownames(retain.cts)<-names(var33.cts);
## setup distribution list for each continuous node
mydists.cts<-list( v2 = "gaussian", v5 = "gaussian",
        v7 = "gaussian", v8 = "gaussian", v13 = "gaussian", 
        v14 = "gaussian", v16 = "gaussian", v17 = "gaussian",
        v22 = "gaussian", v23 = "gaussian", v24 = "gaussian",
        v25 = "gaussian", v29 = "gaussian", v30 = "gaussian",
        v31 = "gaussian", v33 = "gaussian");
## build cache of all the local computations
## this information is needed later when running a model search
mycache.cts<-buildscorecache( data.df=var33.cts,
             data.dists=mydists.cts, dag.banned=banned.cts,
             dag.retained=retain.cts, max.parents=1);


###################################################
### code chunk number 11: abn_v0.85.Rnw:360-366 (eval = FALSE)
###################################################
## # Run a single search heuristic for an additive BN
## heur.res.cts<-search.hillclimber( score.cache=mycache.cts,
##               num.searches=1, seed=0, verbose=FALSE,
##               trace=FALSE, timing.on=FALSE);
## # Setting trace=TRUE, the majority consensus network is
## # plotted as the searches progress


###################################################
### code chunk number 12: abn_v0.85.Rnw:372-403
###################################################
mydag<-matrix( 0, 33, 33); 
colnames(mydag)<-rownames(mydag)<-names(var33);#set names
## create empty DAGs
banned.mix<-matrix( 0, 33, 33);
colnames(banned.mix)<-rownames(banned.mix)<-names(var33);
retain.mix<-matrix( 0, 33, 33);
colnames(retain.mix)<-rownames(retain.mix)<-names(var33);
## setup distribution list for mixed node
mydists.mix<-list( v1 = "binomial", v2 = "gaussian",
        v3 = "binomial", v4 = "binomial", v5 = "gaussian",
        v6 = "binomial", v7 = "gaussian", v8 = "gaussian",
        v9 = "binomial", v10 = "binomial", v11 = "binomial",
        v12 = "binomial", v13 = "gaussian", v14 = "gaussian",
        v15 = "binomial", v16 = "gaussian", v17 = "gaussian",
        v18 = "binomial", v19 = "binomial", v20 = "binomial",
        v21 = "binomial", v22 = "gaussian", v23 = "gaussian",
        v24 = "gaussian", v25 = "gaussian", v26 = "binomial",
        v27 = "binomial", v28 = "binomial", v29 = "gaussian",
        v30 = "gaussian", v31 = "gaussian", v32 = "binomial",
        v33 = "gaussian");
## build cache of all the local computations
## this information is needed later when running a model search
mycache.mix<-buildscorecache( data.df=var33,
             data.dists=mydists.mix, dag.banned=banned.mix,
             dag.retained=retain.mix, max.parents=1);
# Run a single search heuristic for an additive BN
heur.res.mix<-search.hillclimber( score.cache=mycache.mix,
              num.searches=1, seed=0, verbose=FALSE,
              trace=FALSE, timing.on=FALSE);
# Setting trace=TRUE, the majority consensus network is 
# plotted as the searches progress


###################################################
### code chunk number 13: abn_v0.85.Rnw:415-446
###################################################
mydag<-matrix( 0, 33, 33); 
colnames(mydag)<-rownames(mydag)<-names(var33);#set names
## create empty DAGs
banned.mix<-matrix( 0, 33, 33);
colnames(banned.mix)<-rownames(banned.mix)<-names(var33);
retain.mix<-matrix( 0, 33, 33);
colnames(retain.mix)<-rownames(retain.mix)<-names(var33);
## setup distribution list for mixed node
mydists.mix<-list( v1 = "binomial", v2 = "gaussian",
        v3 = "binomial", v4 = "binomial", v5 = "gaussian",
        v6 = "binomial", v7 = "gaussian", v8 = "gaussian",
        v9 = "binomial", v10 = "binomial", v11 = "binomial",
        v12 = "binomial", v13 = "gaussian", v14 = "gaussian",
        v15 = "binomial", v16 = "gaussian", v17 = "gaussian",
        v18 = "binomial", v19 = "binomial", v20 = "binomial",
        v21 = "binomial", v22 = "gaussian", v23 = "gaussian",
        v24 = "gaussian", v25 = "gaussian", v26 = "binomial",
        v27 = "binomial", v28 = "binomial", v29 = "gaussian",
        v30 = "gaussian", v31 = "gaussian", v32 = "binomial",
        v33 = "gaussian");
n.searches<- 10; # example only - must be much larger in practice
## parent limits
max.par<-1 #only 1 because take some minutes for buildscorecache()
## now build cache
mycache.mix<-buildscorecache(data.df=var33, data.dists=mydists.mix,
dag.banned=banned.mix, dag.retained=retain.mix, max.parents=max.par)
# repeat but this time have the majority consensus network plotted 
# as the searches progress
myres.mlp<-search.hillclimber(score.cache=mycache.mix,
           num.searches=n.searches, seed=0, verbose=FALSE,
           trace=FALSE, timing.on=FALSE);


###################################################
### code chunk number 14: abn_v0.85.Rnw:459-464
###################################################
tographviz(dag= myres.mlp$consensus, data.df=var33,
data.dists=mydists.mix, outfile="dagcon.dot");#create file
# dagcon.dot can then be processed with graphviz
# unix shell "dot -Tpdf dagcon.dot -o dagcon.pdf" or use 
# gedit if on Windows


###################################################
### code chunk number 15: abn_v0.85.Rnw:490-531
###################################################
#specific a DAG model - the model we wish to use to perform 
#parametric bootstrapping
mydag.pigs<-matrix(c( 
# D1 D2 D3 D4 D5 D6 D7 D8 D9 D10 Year Loc.x Loc.y
  0, 1, 0, 0, 0, 0, 0, 0, 0, 0,  0,   0,    0,    # D1
  0, 0, 1, 0, 0, 0, 0, 0, 0, 0,  0,   0,    0,    # D2
  0, 0, 0, 1, 0, 0, 0, 0, 0, 0,  0,   0,    0,    # D3
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0,   1,    0,    # D4
  0, 0, 0, 0, 0, 1, 0, 0, 0, 0,  0,   0,    0,    # D5
  0, 0, 0, 1, 0, 0, 0, 0, 0, 0,  0,   0,    0,    # D6
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  1,   0,    0,    # D7
  0, 0, 0, 0, 0, 0, 0, 0, 0, 1,  0,   0,    0,    # D8
  0, 1, 0, 0, 0, 0, 0, 0, 0, 0,  0,   0,    0,    # D9
  0, 0, 0, 0, 0, 0, 0, 0, 1, 0,  0,   0,    0,    # D10
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0,   0,    0,    # Year
  0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0,   0,    0,    # Loc.x
  0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0,   0,    0     # Loc.y
               ),byrow=TRUE,ncol=13); 
colnames(mydag.pigs)<-rownames(mydag.pigs)<-names(pigs.1par);
## setup distribution list for each "pigs" node
mydists.pigs <- list( D1 = "binomial", D2 = "binomial",
        D3 = "binomial", D4 = "binomial", D5 = "binomial",
        D6 = "binomial", D7 = "binomial", D8 = "binomial",
        D9 = "binomial", D10 = "binomial", Year = "gaussian",
        Loc.x = "gaussian", Loc.y = "gaussian");
# Node D1|D2 e.g. logit(P(D1=TRUE)=constant+coeff*D2
# first get the posterior density for the constant and
# get grid of discrete values x and f(x)
marg.D1.1<-fitabn( data.df=pigs.1par, dag.m=mydag.pigs,
           data.dists=mydists.pigs, compute.fixed=TRUE,
           marginal.node=1, marginal.param=1, # intercept
           variate.vec=seq(from=1,to=1.7,len=1000),
           verbose=FALSE, n.grid=1000);
print(names(marg.D1.1$marginals));
# now repeat for the slope term coeff
marg.D1.2<-fitabn( data.df=pigs.1par, dag.m=mydag.pigs,
           data.dists=mydists.pigs, compute.fixed=TRUE,
           marginal.node=1, marginal.param=2, # slope
           variate.vec=seq(from=0.6,to=1.5,len=1000),
           verbose=FALSE, n.grid=1000);
print(names(marg.D1.2$marginals));


###################################################
### code chunk number 16: newfile
###################################################
par( mar=c(8.8,8.2,3.1,3.1),mgp=c(4,2,0));
par( cex.axis=2,cex.lab=2,cex.main=2);
par( las=1,xaxs="i",yaxs="i");
plot( marg.D1.1$marginals[["D1|(Intercept)"]],xlab="Log odds",
      ylab="Density",type="l",axes=!FALSE,xlim=c(0,2), 
      ylim=c(0,7),col="brown",lwd=3,lty=1);
lines( marg.D1.2$marginals[["D1|D2"]],
      col="blue",lwd=3,lty=6);
legend( "topleft",legend=c(expression(paste("f(",beta[D1*","~0],
            "|D)= intercept node D1",sep="")), 
        expression(paste("f(",beta[D1*","~1],
            "|D)= effect of D2 at node D1",sep=""))), cex=1.5, 
        col=c("brown","blue"),lty=c(1,6),lwd=5, bty='n');


###################################################
### code chunk number 17: abn_v0.85.Rnw:582-611
###################################################
pigs<-pigs.1par[,c(1:8,12,13)];
# all 9011 observation but limit to 10 variables
# using all 13 variables in pigs will take several
# hours of cpu time
max.par <- 1
# setup the distribution for each "pigs subset nodes
mydists.pigs.sub <- list( D1 = "binomial", D2="binomial",
        D3 = "binomial", D4 = "binomial", D5 = "binomial",
        D6 = "binomial", D7 = "binomial", D8 = "binomial",
        Loc.x = "gaussian", Loc.y = "gaussian");
banned.pigs.sub<-matrix( 0, 10, 10);#banlist with no constraints
colnames(banned.pigs.sub)<-rownames(banned.pigs.sub)<-names(pigs);
retain.pigs.sub<-matrix( 0, 10, 10);#retainlist without constraints
colnames(retain.pigs.sub)<-rownames(retain.pigs.sub)<-names(pigs);
#compute node cache - note restriction of max. 1 parent 
# per node, this should be increased as necessary
system.time( mynodes.add<-buildscorecache( data.df=pigs,
       data.dists=mydists.pigs.sub, max.parents=max.par,
       dag.banned=banned.pigs.sub, dag.retained=retain.pigs.sub));
# now find the globally best model using previous node cache - so 
# we are only looking for the best DAG, most probable network, 
# within the scope of no more than one parent per node
map.1par.10var<-mostprobable( score.cache=mynodes.add);
tographviz(dag=map.1par.10var,data.df=pigs,
           data.dists=mydists.pigs.sub,
           outfile="map1_10var.dot");#create file
# map1_10var.dot can then be processed with graphviz
# unix shell "dot -Tpdf map1_10var.dot -o map1_10var.pdf" or 
# use gedit if on Windows


###################################################
### code chunk number 18: abn_v0.85.Rnw:622-645
###################################################
pigs.all<-pigs.1par;#all observations all variables
## setup distribution list for each node
mydists.pigs <- list( D1 = "binomial", D2 = "binomial",
        D3 = "binomial", D4 = "binomial", D5 = "binomial",
        D6 = "binomial", D7 = "binomial", D8 = "binomial",
        D9 = "binomial", D10 = "binomial", Year = "gaussian",
        Loc.x = "gaussian", Loc.y = "gaussian");
banned.pigs <- matrix( 0, 13, 13);
colnames(banned.pigs)<-rownames(banned.pigs)<-names(pigs.all);
retain.pigs <-  matrix( 0, 13, 13);
colnames(retain.pigs)<-rownames(retain.pigs)<-names(pigs.all);
system.time( mynodes.add.all<-buildscorecache( data.df=pigs.all,
        max.parents=1, data.dists=mydists.pigs,
        dag.banned=banned.pigs, dag.retained=retain.pigs));
## now for most probable network of all DAGs where each node has 
## at most one arc
system.time( map.1par<-mostprobable( score.cache=mynodes.add.all, 
                                    prior.choice=1));
tographviz( dag=map.1par,data.df=pigs.all, data.dists=mydists.pigs,
           outfile="map_1par.dot");#create file
# mydag.dot can then be processed with graphviz
# unix shell "dot -Tpdf map_1par.dot -o map_1par.pdf" or use gedit 
# if on Windows 


