###############################################################################
## test.R --- 
## Author          : Gilles Kratzer
## Document created: 12/10/2016
## Last modified   : 12/10/2016
###############################################################################
##Purpose: Test the ABN sofware to be inline with previous version + test new updates

Sys.setenv("R_TESTS" = "")
library(abn)
#test_check("abn")
library(testthat)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##Historical tests
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

context("Historical tests")

test_that("Test fitabn()", {
  
  load(file = "testdata/fitabn_ex0.Rdata")
  
  myres.c.test<-fitabn(dag.m = mydag,data.df=mydat,data.dists=mydists,create.graph = T)
  expect_that(myres.c.test,equals(myres.c))
  })

  test_that("Test buildscorecache() and hillclimber()", {
    
    load(file = "testdata/buildscorecache_ex1.Rdata")
    invisible(mycache.test<-buildscorecache(data.df=mydat,data.dists=mydists, max.parents=max.par))
    invisible(mp.dag.test<-mostprobable(score.cache=mycache))
    invisible(myres.test<-fitabn(dag.m=mp.dag,data.df=mydat,data.dists=mydists,create.graph=TRUE))
    invisible(heur.res.test<-search.hillclimber(score.cache=mycache,num.searches=10,seed=42,verbose=F, trace=FALSE,timing.on=TRUE))
  
    expect_that(mycache.test,equals(mycache))
    expect_that(mp.dag.test,equals(mp.dag))
    expect_that(myres.test,equals(myres))
    expect_that(heur.res.test,equals(heur.res))
    
  })
  
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##Updates tests
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  context("Specific tests")


