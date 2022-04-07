###############################################################################
## Purged out of test-package.R
## to have a single context per file as per https://r-pkgs.org/tests.html


context("Historical tests")
Sys.setenv(R_TESTS="")   # really required here ?

test_that("Test fitAbn()", {

    load(file="testdata/fitabn_ex0.Rdata")
    # load(file='tests/testthat/testdata/fitabn_ex0.Rdata')

    myres.c.test <- fitAbn(dag=mydag, data.df=mydat, data.dists=mydists, create.graph=TRUE)

    expect_that(unclass(myres.c.test)[[1]], equals(myres.c[[1]]))
    expect_that(unclass(myres.c.test)[[2]], equals(myres.c[[2]]))
    expect_that(unclass(myres.c.test)[[3]], equals(myres.c[[3]]))
    expect_that(unclass(myres.c.test)[[4]], equals(myres.c[[4]]))
    expect_that(unclass(myres.c.test)[[5]], equals(myres.c[[5]]))
    expect_that(unclass(myres.c.test)[[6]], equals(myres.c[[6]]))
    expect_that(unclass(myres.c.test)[[7]], equals(myres.c[[7]]))
    expect_that(unclass(myres.c.test)[[8]], equals(myres.c[[8]]))
    expect_that(unclass(myres.c.test)[[9]], equals(myres.c[[9]]))
    expect_that(unclass(myres.c.test)[[10]], equals(myres.c[[10]]))
    expect_that(unclass(myres.c.test)[[11]], equals(myres.c[[11]]))
    expect_that(unclass(myres.c.test)[[12]], equals(myres.c[[12]]))
    expect_that(unclass(myres.c.test)[[13]], equals(myres.c[[13]]))
#     expect_that(unclass(myres.c.test)[[14]], equals(myres.c[[14]])) # no graph  part anymore...

})

test_that("Test buildScoreCache.bayes() and hillclimber()", {
#    if(requireNamespace("INLA", quietly=TRUE)){
    load(file="testdata/buildscorecache_ex1.Rdata")
    # load(file='tests/testthat/testdata/buildscorecache_ex1.Rdata')

#    invisible(mycache.test <- abn:::buildScoreCache.bayes(data.df=mydat, data.dists=mydists, max.parents=max.par)) ??
    invisible(mycache.test <- abn:::buildScoreCache.bayes(data.df=mydat, data.dists=mydists, max.parents=max.par))
    class(mycache.test) <-  c("abnCache")
    invisible(mp.dag.test <- mostProbable(score.cache=mycache.test, verbose=FALSE))
    invisible(myres.test <- fitAbn(dag=mp.dag, data.df=as.data.frame(mydat), data.dists=mydists, create.graph=TRUE))
    invisible(heur.res.test <- searchHillClimber(score.cache=mycache.test, num.searches=10, seed=42, verbose=FALSE, timing.on=TRUE))


    expect_that(mycache.test[1:3], equals(mycache[1:3]))
    expect_that(unclass(mp.dag.test[[1]]), equals((mp.dag)))

    expect_that(unclass(myres.test)[[1]], equals(myres[[1]]))
    expect_that(unclass(myres.test)[[2]], equals(myres[[2]]))
    expect_that(unclass(myres.test)[[3]], equals(myres[[3]]))
    expect_that(unclass(myres.test)[[4]], equals(myres[[4]]))
    expect_that(unclass(myres.test)[[5]], equals(myres[[5]]))
    expect_that(unclass(myres.test)[[6]], equals(myres[[6]]))
    expect_that(unclass(myres.test)[[7]], equals(myres[[7]]))
    expect_that(unclass(myres.test)[[8]], equals(myres[[8]]))
    expect_that(unclass(myres.test)[[9]], equals(myres[[9]]))
    expect_that(unclass(myres.test)[[10]], equals(myres[[10]]))
    expect_that(unclass(myres.test)[[11]], equals(myres[[11]]))
    expect_that(unclass(myres.test)[[12]], equals(myres[[12]]))
    expect_that(unclass(myres.test)[[13]], equals(myres[[13]]))
    expect_that(unclass(myres.test)[[14]], equals(myres[[14]]))

    expect_that(heur.res.test[[1]], equals(heur.res[[1]]))
    expect_that(heur.res.test[[2]], equals(heur.res[[2]]))
    expect_that(heur.res.test[[3]], equals(heur.res[[3]]))
    expect_that(heur.res.test[[4]], equals(heur.res[[4]]))
    expect_that(heur.res.test[[5]], equals((heur.res[[5]])))
    expect_that(heur.res.test[[6]], equals(heur.res[[6]]))
#    } else {
#       cat('Package INLA not available, number of passed tests might be different')
#    }
})

