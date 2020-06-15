############################################################################### 
## test.R --- Author : Gilles Kratzer 
## Document created: 29/05/2019 
## Last modified :  
##  Purpose: Test the ABN sofware to be maintain backward compatibility

#Sys.setenv(R_TESTS = "")

context("loading library")

# test_check('abn')

#suppressWarnings(library(abn))

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
## Historical tests ()
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

context("Test backward compatibility")

test_that("buildscorecache", {
  
  test_unit <- function(){
    
    mydat <- ex0.dag.data[,c("b1","b2","g1","g2","b3","g3")];## take a subset of cols
    
    ## setup distribution list for each node
    mydists <- list(b1="binomial",
                  b2="binomial",
                  g1="gaussian",
                  g2="gaussian",
                  b3="binomial",
                  g3="gaussian"
    );
    
    ban <- matrix(rep(0,dim(mydat)[2]^2),ncol=dim(mydat)[2]);# ban nothing
    colnames(ban) <- rownames(ban) <- names(mydat); #names must be set
    ban["b1","b2"] <- 1; # now ban arc from b2 to b1 
    retain <- matrix(rep(0,dim(mydat)[2]^2),ncol=dim(mydat)[2]);# retain nothing
    colnames(retain) <- rownames(retain) <- names(mydat); #names must be set
    retain["g1","g3"] <- 1; # always retain arc from g3 to g1
    # parent limits
    max.par <- list("b1"=2,"b2"=2,"g1"=2,"g2"=0,"b3"=2,"g3"=3);
    
    ## now build cache of scores (goodness of fits for each node)
    
    res.c <- buildscorecache(data.df=mydat,data.dists=mydists,
                           dag.banned=ban, dag.retained=retain,max.parents=max.par
    );
    
    ## repeat but using R-INLA. The mlik's should be virtually identical.
    ## now build cache
    res.inla <- buildscorecache(data.df=mydat,data.dists=mydists,
                              dag.banned=ban, dag.retained=retain,max.parents=max.par,
                              max.mode.error=100);
    
    
    
    #################################################################
    ## Example 2 - much bigger problem using glm - make take a while
    #################################################################
    
    mydat <- ex2.dag.data;## this data comes with abn see ?ex2.dag.data
    
    ## setup distribution list for each node
    mydists <- list(b1="binomial",
                  g1="gaussian",
                  p1="poisson",
                  b2="binomial",
                  g2="gaussian",
                  p2="poisson",
                  b3="binomial",
                  g3="gaussian",
                  p3="poisson",
                  b4="binomial",
                  g4="gaussian",
                  p4="poisson",
                  b5="binomial",
                  g5="gaussian",
                  p5="poisson",
                  b6="binomial",
                  g6="gaussian",
                  p6="poisson"
    );
    
    ## parent limits
    max.par <- list("b1"=2,"g1"=2,"p1"=2,"b2"=2,"g2"=2,"p2"=2,"b3"=2,
                  "g3"=2,"p3"=2,"b4"=2,"g4"=2,
                  "p4"=2,"b5"=2,"g5"=2,"p5"=2,"b6"=2,"g6"=2,"p6"=2);
    
    ## no explicit ban or retain restrictions set so dont need to supply ban
    ##  or retain matrices
    
    ## now build cache using internal code just for nodes 1,2 and 3
    ## e.g. "b1", "p1" and "g1" 
    mycache.c <- buildscorecache(data.df=mydat,data.dists=mydists,
                               max.parents=max.par, which.nodes=c(1:3));
    
    ###################################################################
    ## Example 3 - grouped data - random effects example e.g. glmm
    ###################################################################
    
    mydat <- ex3.dag.data;## this data comes with abn see ?ex3.dag.data
    
    mydists <- list(b1="binomial",
                  b2="binomial",
                  b3="binomial",
                  b4="binomial",
                  b5="binomial",
                  b6="binomial",
                  b7="binomial",
                  b8="binomial",
                  b9="binomial",
                  b10="binomial",
                  b11="binomial",
                  b12="binomial",
                  b13="binomial"
    );
    max.par <- 1;
    
    ## in this example INLA is used as default since these are glmm nodes
    ## when running this at node-parent combination 71 the default accuracy check on the 
    ## INLA modes is exceeded (default is a max. of 10 percent difference from
    ## modes estimated using internal code) and a message is given that internal code
    ## will be used in place of INLA's results.
    
    mycache <- buildscorecache(data.df=mydat,data.dists=mydists,group.var="group",
                             cor.vars=c("b1","b2"),
                             max.parents=max.par, which.nodes=c(1));
    
    mydat <- ex0.dag.data[,c("b1","b2","g1","g2","b3","g3")] ## take a subset of cols
    
    ## setup distribution list for each node
    mydists <- list(b1="binomial",
                    b2="binomial",
                    g1="gaussian",
                    g2="gaussian",
                    b3="binomial",
                    g3="gaussian")
    
    ## now build cache of scores (goodness of fits for each node)
    res.mle <- buildscorecache(data.df=mydat,data.dists=mydists,max.parents=3,method="mle")
    res.abn <- buildscorecache(data.df=mydat,data.dists=mydists,max.parents=3,method="Bayes")
    
  }
  if(requireNamespace("INLA", quietly = TRUE)){
  expect_error(test_unit(),NA)
  }
})


test_that("mostprobable", {
  
  test_unit <- function(){
    
    ## use built-in simulated data set

    mydat <- ex0.dag.data[,c("b1","b2","g1","g2","p1","p2")];
    ## take a subset of cols

    ## setup distribution list for each node
    mydists <- list(b1="binomial",
                  b2="binomial",
                  g1="gaussian",
                  g2="gaussian",
                  p1="poisson",
                  p2="poisson"
    );

    #use simple banlist with no constraints
    ban <- matrix(c(
      #   1 2 3 4 5 6
      0,0,0,0,0,0, # 1
      0,0,0,0,0,0, # 2
      0,0,0,0,0,0, # 3
      0,0,0,0,0,0, # 4
      0,0,0,0,0,0, # 5
      0,0,0,0,0,0 # 6
    ),byrow=TRUE,ncol=6);

    colnames(ban) <- rownames(ban) <- names(mydat); #names must be set
    ban["b1","b2"] <- 1; # now ban arc from b2 to b1

    retain <- matrix(c(
      #   1 2 3 4 5 6
      0,0,0,0,0,0, # 1
      0,0,0,0,0,0, # 2
      0,0,0,0,0,0, # 3
      0,0,0,0,0,0, # 4
      0,0,0,0,0,0, # 5
      0,0,0,0,0,0 # 6
    ),byrow=TRUE,ncol=6);

    colnames(retain) <- rownames(retain) <- names(mydat); #names must be set
    retain["g1","g2"] <- 1; # always retain arc from g2 to g1
    # parent limits
    max.par <- list("b1"=1,"b2"=1,"g1"=1,"g2"=0,"p1"=1,"p2"=2);
    ## now build cache
    mycache <- buildscorecache(data.df=mydat,data.dists=mydists,
                             dag.banned=ban, dag.retained=retain,max.parents=max.par);

    #now find the globally best DAG
    mp.dag <- mostprobable(score.cache=mycache, verbose = FALSE);
    # get the corresponding best goodness of fit - network score
    #m <- fitabn(dag.m=mp.dag,data.df=mydat,data.dists=mydists)$mlik;
    
    ## Second example ############
    
    mydat <- ex1.dag.data;## this data comes with abn see ?ex1.dag.data
    
    ## setup distribution list for each node
    mydists <- list(b1="binomial",
                  p1="poisson",
                  g1="gaussian",
                  b2="binomial",
                  p2="poisson",
                  b3="binomial",
                  g2="gaussian",
                  b4="binomial",
                  b5="binomial",
                  g3="gaussian"
    );
    
    #use simple banlist with no constraints
    ban <- matrix(c(
      #   1 2 3 4 5 6    
      0,0,0,0,0,0,0,0,0,0, # 1 
      0,0,0,0,0,0,0,0,0,0, # 2
      0,0,0,0,0,0,0,0,0,0, # 3 
      0,0,0,0,0,0,0,0,0,0, # 4
      0,0,0,0,0,0,0,0,0,0, # 5 
      0,0,0,0,0,0,0,0,0,0, # 6
      0,0,0,0,0,0,0,0,0,0, # 7 
      0,0,0,0,0,0,0,0,0,0, # 8
      0,0,0,0,0,0,0,0,0,0, # 9 
      0,0,0,0,0,0,0,0,0,0  # 10     
    ),byrow=TRUE,ncol=10);
    
    colnames(ban) <- rownames(ban) <- names(mydat); #names must be set
    
    retain <- matrix(c(
      #   1 2 3 4 5 6    
      0,0,0,0,0,0,0,0,0,0, # 1 
      0,0,0,0,0,0,0,0,0,0, # 2
      0,0,0,0,0,0,0,0,0,0, # 3 
      0,0,0,0,0,0,0,0,0,0, # 4
      0,0,0,0,0,0,0,0,0,0, # 5 
      0,0,0,0,0,0,0,0,0,0, # 6
      0,0,0,0,0,0,0,0,0,0, # 7 
      0,0,0,0,0,0,0,0,0,0, # 8
      0,0,0,0,0,0,0,0,0,0, # 9 
      0,0,0,0,0,0,0,0,0,0  # 10     
    ),byrow=TRUE,ncol=10);
    colnames(retain) <- rownames(retain) <- names(mydat); #names must be set
    
    ## parent limits
    max.par <- list("b1"=2,"p1"=2,"g1"=2,"b2"=2,"p2"=2,"b3"=2,"g2"=2,"b4"=2,"b5"=2,"g3"=2);
    ## now build cache
    mycache.exemple1 <- buildscorecache(data.df=mydat,data.dists=mydists,
                             dag.banned=ban, dag.retained=retain,max.parents=max.par);
    
    #now find the globally best DAG
    mp.dag <- mostprobable(score.cache=mycache.exemple1, verbose = FALSE);
    #m <- fitabn(dag.m=mp.dag,data.df=mydat,data.dists=mydists)$mlik;
    
    ## plot the best model
    #myres <- fitabn(dag.m=mp.dag,data.df=mydat,data.dists=mydists,create.graph=TRUE);
    
    
    ## fit the known true DAG
    # true.dag <- ban;
    # true.dag["p2",c("b1","p1")] <- 1;
    # true.dag["b3",c("b1","g1","b2")] <- 1;
    # true.dag["g2",c("p1","g1","b2")] <- 1;
    # true.dag["b4",c("g1","p2")] <- 1;
    # true.dag["b5",c("g1","g2")] <- 1;
    # true.dag["g3",c("g1","b2")] <- 1;
    # 
    # m <- fitabn(dag.m=true.dag,data.df=mydat,data.dists=mydists)$mlik;
    
    #################################################################
    ## example 3 - models with random effects
    #################################################################
    
    mydat <- ex3.dag.data[,c(1:4,14)];
    ## this data comes with abn see ?ex3.dag.data
    
    mydists <- list(b1="binomial",
                  b2="binomial",
                  b3="binomial",
                  b4="binomial"
    );
    max.par <- 1;
    
    mycache.mixed <- buildscorecache(data.df=mydat,data.dists=mydists,
                                   group.var="group",cor.vars=c("b1","b2","b3","b4"),
                                   ## each node uses a random effect adjustment
                                   max.parents=max.par);                         
    
    ## find the most probable DAG
    #mp.dag <- mostprobable(score.cache=mycache.c);
    
    ## get goodness of fit
    #m <- fitabn(dag.m=mp.dag,data.df=mydat,data.dists=mydists,group.var="group",cor.vars=c("b1","b2","b3","b4"))$mlik;
  }
  if(requireNamespace("INLA", quietly = TRUE)){
  expect_error(test_unit(),NA)
  }
})

test_that("search.hillclimber", {
  
  test_unit <- function(){
    ##############################################
    ## example 1: use built-in simulated data set
    ##############################################
    
    mydat <- ex1.dag.data; ## this data comes with abn see ?ex1.dag.data

    ## setup distribution list for each node
    mydists <- list(b1="binomial",
                  p1="poisson",
                  g1="gaussian",
                  b2="binomial",
                  p2="poisson",
                  b3="binomial",
                  g2="gaussian",
                  b4="binomial",
                  b5="binomial",
                  g3="gaussian"
    );

    ## not run because may take some minutes for buildscorecache()
    ## parent limits
    max.par <- list("b1"=2,"p1"=2,"g1"=2,"b2"=2,"p2"=2,"b3"=2,"g2"=2,"b4"=2,"b5"=2,"g3"=2);
    ## now build cache

    mycache <- buildscorecache(data.df=mydat,data.dists=mydists,max.parents=max.par);

    heur.res <- search.hillclimber(score.cache=mycache,
                                 num.searches=100,seed=0,verbose=FALSE,timing.on=FALSE);
    
    ###################################################################################################
    ## example 2 - glmm example - but no difference here as the format of the score cache is identical
    ###################################################################################################
    
    mydat <- ex3.dag.data[,c(1:5,14)];## this data comes with abn see ?ex1.dag.data

    mydists <- list(b1="binomial",
                  b2="binomial",
                  b3="binomial",
                  b4="binomial",
                  b5="binomial"
    );
    max.par <- 1;

    mycache.mixed <- buildscorecache(data.df=mydat,data.dists=mydists,group.var="group",
                                   cor.vars=c("b1","b2","b3","b4","b5"),
                                   max.parents=max.par, which.nodes=c(1:5));
    
    # now peform 1000 greedy searches
    heur.res <- search.hillclimber(score.cache=mycache.mixed,num.searches=100,
                                 seed=0,verbose=FALSE,timing.on=FALSE);
    
  }
  if(requireNamespace("INLA", quietly = TRUE)){
  expect_error(test_unit(),NA)
  }
})


test_that("fitabn", {

    test_unit <- function(){
        ## use built-in simulated data set
        
        mydat <- ex0.dag.data[,c("b1","b2","b3","g1","b4","p2","p4")];## take a subset of cols
        
        ## setup distribution list for each node
        mydists <- list(b1="binomial",
                      b2="binomial",
                      b3="binomial",
                      g1="gaussian",
                      b4="binomial",
                      p2="poisson",
                      p4="poisson"
        );
        ## null model - all independent variables
        mydag.empty <- matrix(data=c(
            0,0,0,0,0,0,0, # 
            0,0,0,0,0,0,0, #
            0,0,0,0,0,0,0, # 
            0,0,0,0,0,0,0, # 
            0,0,0,0,0,0,0, #
            0,0,0,0,0,0,0, #
            0,0,0,0,0,0,0  #
        ), byrow=TRUE,ncol=7);
        colnames(mydag.empty) <- rownames(mydag.empty) <- names(mydat);
        
        ## now fit the model to calculate its goodness of fit
        myres.c <- fitabn(dag.m=mydag.empty,data.df=mydat,data.dists=mydists,centre=TRUE,
                        compute.fixed=FALSE);
        
         
        ## now repeat but include some dependencies
        mydag <- mydag.empty;
        mydag["b1","b2"] <- 1; # b1 <- b2 
        mydag["b2","p4"] <- 1; # b2 <- p4
        mydag["b2","g1"] <- 1; # b2 <- g1
        mydag["g1","p2"] <- 1; # g1 <- p2
        mydag["b3","g1"] <- 1; # b3 <- g1
        mydag["b4","b1"] <- 1; # b4 <- b1
        mydag["p4","g1"] <- 1; # p4 <- g1
        
        myres.c <- fitabn(dag.m=mydag,data.df=mydat,data.dists=mydists,centre=TRUE,
                        compute.fixed=FALSE);
        
       
        ## now also plot the graph via graphviz 
        myres.c <- fitabn(dag.m=mydag,data.df=mydat,data.dists=mydists,centre=TRUE,
                        create.graph=TRUE,compute.fixed=FALSE);
        
       
        ## a simple plot of some posterior densities
        ## the algorithm which chooses density points is very simple any may be
        ## rather sparse so also recompute the density over an equally spaced
        ## grid of 100 points between the two end points
        ## which had at f=min.pdf
        myres.c <- fitabn(dag.m=mydag,data.df=mydat,data.dists=mydists,centre=TRUE,
                        compute.fixed=TRUE,n.grid=100);
        
        ## repeat but use INLA for the numerics using max.mode.error=100
        ## as using internal code is the default here rather than INLA
        myres.inla <- fitabn(dag.m=mydag,data.df=mydat,data.dists=mydists,centre=TRUE,
                           compute.fixed=TRUE,max.mode.error=100);
        
      ### A very simple mixed model example using built-in data
        ## specify DAG - only two variables using subset of variables from ex3.dag.data
        ## both variables are assumed to need (separate) adjustment for the group variable
        ## i.e. a binomial glmm at each node
        
        ## model where b1 <- b2
        mydag <- matrix(data=c(
            0,1, # b1
            0,0  # b2
        ), byrow=TRUE,ncol=2);
        
        colnames(mydag) <- rownames(mydag) <- names(ex3.dag.data[,c(1,2)]);
        ## specific distributions
        mydists <- list(b1="binomial",
                      b2="binomial"
        );
        
        ## just compute marginal likelihood - use internal code via max.mode.error=0
        ## as using INLA is the default here.
        myres.c <- fitabn(dag.m=mydag,data.df=ex3.dag.data[,c(1,2,14)],data.dists=mydists,
                        group.var="group",cor.vars=c("b1","b2"),
                        centre=TRUE,compute.fixed=FALSE,max.mode.error=0);
         
        ## compare with INLA estimate
        myres.inla <- fitabn(dag.m=mydag,data.df=ex3.dag.data[,c(1,2,14)],
                           data.dists=mydists,group.var="group",cor.vars=c("b1","b2"),
                           centre=TRUE,compute.fixed=FALSE,max.mode.error=100); 
        
        ## compare log marginal likelihoods for each node and total DAG - should be very similar
        
        ## now for marginals - INLA is strongly preferable for estimating marginals for nodes 
        ## with random effects as it is far faster, but may not be reliable
        ## see www.r-bayesian-networks.org/quality-assurance
        
        # INLA's estimates of the marginals, using high n.grid=500
        # as this makes the plots smoother - see below.
        # myres.inla <- fitabn(dag.m=mydag,data.df=ex3.dag.data[,c(1,2,14)],data.dists=mydists,
        #                    group.var="group",cor.vars=c("b1","b2"),
        #                    compute.fixed=TRUE,max.mode.error=100,
        #                    n.grid=500,max.hessian.error = 10E-02);
        
        ## this is NOT recommended - marginal density estimation using fitabn in mixed models
        ## is really just for diagnostic purposes, better to use fitabn.inla() here
        ## but here goes...be patient
        # myres.c <- fitabn(dag.m=mydag,data.df=ex3.dag.data[,c(1,2,14)],data.dists=mydists,
        #                 group.var="group",cor.vars=c("b1","b2"),compute.fixed=TRUE,
        #                 max.mode.error=0,max.hessian.error = 10E-02);
        
         
        ### these are very similar although not exactly identical
        
        ## use internal code but only to compute a single parameter over a specified grid
        ## this can be necessary if the simple auto grid finding functions does a poor job
        
        # myres.c <- fitabn(dag.m=mydag,data.df=ex3.dag.data[,c(1,2,14)],data.dists=mydists,
        #                 group.var="group",
        #                 cor.vars=c("b1","b2"),centre=FALSE,compute.fixed=TRUE,
        #                 marginal.node=1,marginal.param=3,## precision term in node 1
        #                 variate.vec=seq(0.05,1.5,len=25),max.hessian.error = 10E-02);
        
  
    ## use built-in simulated data set
    
    mydat <- ex0.dag.data[,c("b1","b2","b3","g1","b4","p2","p4")];## take a subset of cols
    
    ## setup distribution list for each node
    mydists <- list(b1="binomial",
                    b2="binomial",
                    b3="binomial",
                    g1="gaussian",
                    b4="binomial",
                    p2="poisson",
                    p4="poisson"
    );
    ## null model - all independent variables
    mydag.empty <- matrix(data=c(
        0,0,0,0,0,0,0, # 
        0,0,0,0,0,0,0, #
        0,0,0,0,0,0,0, # 
        0,0,0,0,0,0,0, # 
        0,0,0,0,0,0,0, #
        0,0,0,0,0,0,0, #
        0,0,0,0,0,0,0  #
    ), byrow=TRUE,ncol=7);
    colnames(mydag.empty) <- rownames(mydag.empty) <- names(mydat);
    ## now repeat but include some dependencies
    mydag <- mydag.empty;
    mydag["b1","b2"] <- 1; # b1 <- b2 
    mydag["b2","p4"] <- 1; # b2 <- p4
    mydag["b2","g1"] <- 1; # b2 <- g1
    mydag["g1","p2"] <- 1; # g1 <- p2
    mydag["b3","g1"] <- 1; # b3 <- g1
    mydag["b4","b1"] <- 1; # b4 <- b1
    mydag["p4","g1"] <- 1; # p4 <- g1
    
    ## now fit the model to calculate its goodness of fit
    myres.c <- fitabn(dag.m=mydag,data.df=mydat,data.dists=mydists,method="mle",centre=TRUE)
    
    }
    if(requireNamespace("INLA", quietly = TRUE)){
    expect_error(test_unit(),NA)
    }
})


