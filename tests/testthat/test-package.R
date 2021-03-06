###############################################################################
## test.R --- Author : Gilles Kratzer
## Document created: 12/10/2016
## Last modified : 12/10/2016 :
##                 06/12/2016 :
##                 29/04/2017 (update tests):
##                 29/08/2017 (update tests: skweness, fitabn.mle, buildScoreCache.mle, Markov Blanket, Formula statement,infoDag(),essentialgraph(),or(), mi.data(),discretization(), entropy.data(),logit(), expit()) :
##                 16/07/2018(update tests: compareDag. Renaming) :
##                 21/05/2019 (gk: update buildScoreCache tests)
##  Purpose: Test the ABN sofware to be  inline with previous version + test new updates

context("Development tests")

test_that("Test mostProbable()", {

    load(file="testdata/buildscorecache_ex1.Rdata")
    # load(file='tests/testthat/testdata/buildscorecache_ex1.Rdata')
#    if(requireNamespace("INLA", quietly=TRUE)){
      invisible(mycache.test <- abn:::buildScoreCache.bayes(data.df=mydat, data.dists=mydists, max.parents=1))
      class(mycache.test) <-  c("abnCache")
      expect_silent(mp.dag.test <- mostProbable(score.cache=mycache.test, verbose=FALSE))
#    } else {
#        cat('Package INLA not available, number of passed tests might be different')
#        expect_equal(1,1)
#    }
})

test_that("Test hill.climber()", {

    load(file="testdata/buildscorecache_ex1.Rdata")
    # load(file='tests/testthat/testdata/buildscorecache_ex1.Rdata')
#    if(requireNamespace("INLA", quietly=TRUE)){
      invisible(mycache.test <- abn:::buildScoreCache.bayes(data.df=mydat, data.dists=mydists, max.parents=1))
      class(mycache.test) <-  c("abnCache")
      class(mycache) <-  c("abnCache")      # why is mychache.test constructed here??
      expect_silent(mp.dag.test <- searchHillClimber(score.cache=mycache, verbose=FALSE))
#    }
#    expect_equal(1,1)
})

test_that("Test fitAbn.mle()", {

    # Gaussian
    df <- airquality[complete.cases(airquality), ]

    dist <- list(Ozone="gaussian", Solar.R="gaussian", Wind="gaussian", Temp="gaussian", Month="gaussian", Day="gaussian")
    names(dist) <- colnames(df)

    d <- matrix(data=0, nrow=6, ncol=6)
    d[1, ] <- c(0, 1, 1, 1, 1, 1)
    colnames(d) <- rownames(d) <- names(dist)

    ## test wrapper

    m.0.mle <- abn:::fitAbn.mle(dag=d, data.df=df, data.dists=dist)
    m.0.mle.1 <- fitAbn(dag=d, data.df=df, data.dists=dist, method="mle")

 #   if(requireNamespace("INLA", quietly=TRUE)){
    m.0.bayes <- abn:::fitAbn.bayes(dag=d, data.df=df, data.dists=dist)
    m.0.bayes.1 <- fitAbn(dag=d, data.df=df, data.dists=dist, method="bayes")
 #   } # else expect_equal(1,1)
    mydists <- list(b1="binomial",  b2="binomial")

    ## test formula statement

 #   expect_error(fitAbn(dag=~b1|b2, data.df=ex3.dag.data[,c(1,2)], data.dists=mydists, max.mode.error=0), NA)
    expect_silent(fitAbn(dag=~b1|b2, data.df=ex3.dag.data[,c(1,2)], data.dists=mydists, max.mode.error=0))

    if(requireNamespace("INLA", quietly=TRUE)){
    ## test formula statement in bayes with correlation
      expect_silent(fitAbn(dag=~b1|b2, data.df=ex3.dag.data[,c(1,2,14)], data.dists=mydists,
                      group.var="group", cor.vars=c("b1","b2"), max.mode.error=0))
    } else expect_equal(1,1)
    expect_silent(fitAbn(dag=~b1|b2, data.df=ex3.dag.data[,c(1,2)], data.dists=mydists, method="mle"))


    ## test wrapper fitAbn method
    expect_equal(m.0.mle,(m.0.mle.1))
#    if(requireNamespace("INLA", quietly=TRUE)){
      expect_equal(unclass(m.0.bayes),unclass(m.0.bayes.1))
#    } #else expect_equal(1,1)

    m1 <- abn:::fitAbn.mle(dag=d, data.df=df, data.dists=dist, centre=FALSE)
    m2 <- lm(df[, 1] ~ as.matrix(df[, 2:6]))

    expect_that(unname(m1$coef[[1]]), equals(unname(t(coef(summary.lm(object=m2))[, 1]))))
    expect_that(unname(m1$Stderror[[1]]), equals(unname(t(coef(summary.lm(object=m2))[, 2]))))
    expect_that(unname(m1$pvalue[[1]]), equals(unname(t(coef(summary.lm(object=m2))[, 4]))))

    ## test centre
    m1 <- abn:::fitAbn.mle(dag=d, data.df=df, data.dists=dist, centre=TRUE)
    m3 <- abn:::fitAbn.mle(dag=d, data.df=df, data.dists=dist)
    d[1, ] <- c(0, 1, 0, 0, 0, 0)
    m2 <- abn:::fitAbn.mle(dag=d, data.df=df, data.dists=dist)
    m4 <- cor(df[, 1:6])

    expect_that(m1, equals(m3))
    expect_that(unname(m2$coef[[1]])[2], equals(m4[1, 2]))

    # Binomial
    dist <- list(a="binomial", b="binomial")

    data.param <- matrix(data=c(0, 0.5, 0, 0), nrow=2L, ncol=2L, byrow=TRUE)

    # naming
    colnames(data.param) <- rownames(data.param) <- names(dist)

    out.sim <- (simulateAbn(data.dists=dist, n.chains=1, n.adapt=100, n.thin=1, n.iter=100, data.param=data.param,
        simulate=TRUE, seed=132,verbose=FALSE))

    m1 <- abn:::fitAbn.mle(dag=data.param, data.df=out.sim, data.dists=dist, centre=FALSE)
    m2 <- glm(formula=out.sim$a ~ out.sim$b, family="binomial")

    expect_that(unname(m1$coef[[1]]), equals(unname(t(coef(summary.glm(object=m2))[, 1]))))
    expect_that(unname(m1$Stderror[[1]]), equals(unname(t(coef(summary.glm(object=m2))[, 2]))))
    expect_that(unname(m1$pvalue[[1]]), equals(unname(t(coef(summary.glm(object=m2))[, 4]))))

    # Poisson
    dist <- list(a="poisson", b="poisson")

    data.param <- matrix(data=c(0, 0.5, 0, 0), nrow=2L, ncol=2L, byrow=TRUE)

    # naming
    colnames(data.param) <- rownames(data.param) <- names(dist)

    out.sim <- invisible(simulateAbn(data.dists=dist, n.chains=1, n.adapt=100, n.thin=1, n.iter=100, data.param=data.param,
        simulate=TRUE, seed=132,verbose=FALSE))

    m1 <- abn:::fitAbn.mle(dag=data.param, data.df=out.sim, data.dists=dist, centre=FALSE)
    m2 <- glm(formula=out.sim$a ~ out.sim$b, family="poisson")

    ## pvalues and stderr are computed up to 10e-06 precision!
    expect_that(unname(m1$coef[[1]]), equals(unname(t(coef(summary.glm(object=m2))[, 1]))))
    expect_that(object=unname(m1$Stderror[[1]]), condition=equals(unname(t(coef(summary.glm(object=m2))[, 2])), tolerance=1e-06))
    expect_that(unname(m1$pvalue[[1]]), equals(unname(t(coef(summary.glm(object=m2))[, 4])), tolerance=1e-06))

    # multinomial (as response and as response)
    dist <- list(a="multinomial", b="gaussian")

    # link matrix
    data.param <- matrix(data=c(1, 0, 0, 1), nrow=2L, ncol=2L, byrow=TRUE)
    # naming
    colnames(data.param) <- rownames(data.param) <- names(dist)

    data.param.mult <- matrix(data=c(0.7, 0.1, 0.2, 0, 0, 0), nrow=2L, ncol=3L, byrow=TRUE)

    out <- invisible(simulateAbn(data.dists=dist, n.chains=1, n.adapt=100, n.thin=1, n.iter=1000, data.param=data.param,
        simulate=TRUE, data.param.mult=data.param.mult, seed=132,verbose=FALSE))

    res <- out$`a[1]`
    res[out$`a[2]` == 1] <- 2
    res[out$`a[3]` == 1] <- 3

    dta <- data.frame(b=rnorm(n=length(res), mean=1, sd=1), c=rnorm(n=length(res), mean=5, sd=1))

    dta$res <- factor(res)

    dist <- list(a="gaussian", b="gaussian", c="multinomial")

    dag <- matrix(data=c(0, 1, 1, 0, 0, 0, 0, 0, 0), nrow=3, ncol=3, byrow=TRUE)

    colnames(dag) <- rownames(dag) <- names(dist)

    names(dta) <- names(dist)

    fit.abn.mle <- abn:::fitAbn.mle(dag=(dag), data.df=dta, data.dists=dist, centre=FALSE)

    m.1 <- model.matrix(object=~dta$c + 0)#, contrasts.arg=list('dta$c'="contr.treatment"))

    m2 <- lm(dta[, 1] ~ -1 + as.matrix(cbind(dta[, 2], m.1)))

    # as predictors

    expect_that(unname(fit.abn.mle$coef[[1]]), equals(unname(t(coef(summary.lm(object=m2))[, 1]))))
    expect_that(unname(fit.abn.mle$Stderror[[1]]), equals(unname(t(coef(summary.lm(object=m2))[, 2])), tolerance=0.001))
    expect_that(unname(fit.abn.mle$pvalue[[1]]), equals(unname(t(coef(summary.lm(object=m2))[, 4])), tolerance=0.001))

    # as response (intercept only)
    expect_that(unname(fit.abn.mle$coef[[3]]), equals(unname(t(coef(summary(multinom(formula=m.1 ~ 1, trace=FALSE)))[, 1]))))
    expect_that(unname(fit.abn.mle$Stderror[[3]]), equals(unname(t(summary(multinom(formula=m.1 ~ 1, trace=FALSE))$standard.errors)),
        tolerance=1e-06))

    z <- summary(multinom(formula=m.1 ~ 1, trace=FALSE))$coefficients/summary(multinom(formula=m.1 ~ 1, trace=FALSE))$standard.errors
    p <- (1 - pnorm(abs(z), 0, 1)) * 2

    expect_that(unname(fit.abn.mle$pvalue[[3]]), equals(unname(t(p)),
        tolerance=1e-10))

    # as response

    dag <- matrix(data=c(0, 0, 0, 0, 0, 0, 1, 1, 0), nrow=3, ncol=3, byrow=TRUE)

    colnames(dag) <- rownames(dag) <- names(dist)

    fit.abn.mle <- abn:::fitAbn.mle(dag=(dag), data.df=dta, data.dists=dist, centre=FALSE)

    m.1 <- model.matrix(object=~dta$c + 0)#, contrasts.arg=contr.treatment)

    m2 <- multinom(formula=m.1 ~ dta$a + dta$b, trace=FALSE)

    expect_that(unname(as.vector(fit.abn.mle$coef[[3]])), equals(unname(as.vector(coef(summary(multinom(formula=m.1 ~ dta$a +
        dta$b, trace=FALSE)))))))
    expect_that(unname(as.vector(fit.abn.mle$Stderror[[3]])), equals(unname(as.vector(summary(multinom(formula=m.1 ~ dta$a +
        dta$b, trace=FALSE))$standard.errors)), tolerance=1e-06))

    # 2-tailed Wald z tests to test significance of coefficients
    z <- summary(multinom(formula=m.1 ~ dta$a + dta$b, trace=FALSE))$coefficients/summary(multinom(formula=m.1 ~ dta$a + dta$b, trace=FALSE))$standard.errors
    p <- (1 - pnorm(abs(z), 0, 1)) * 2

    expect_that(unname(as.vector(fit.abn.mle$pvalue[[3]])), equals(unname(as.vector(p)), tolerance=1e-06))

})

test_that("Test buildScoreCache()", {

    df <- airquality[complete.cases(airquality), ]

    # distribution (gaussian)
    dist <- list(Ozone="gaussian", Solar.R="gaussian", Wind="gaussian", Temp="gaussian", Month="gaussian", Day="gaussian")
    names(dist) <- colnames(df)



    mycache.computed.mle <- abn:::buildScoreCache.mle(data.df=df, data.dists=dist, max.parents=6)
    mycache.old <- abn:::buildScoreCache.bayes(data.df=df, data.dists=dist, max.parents=6, dry.run=TRUE)

    mycache.computed.mle.1 <- abn:::buildScoreCache.mle(data.df=df, data.dists=dist, max.parents=3)
    mycache.old.1 <- abn:::buildScoreCache.bayes(data.df=df, data.dists=dist, max.parents=3, dry.run=TRUE)
    mycache.computed.mle.1.1 <- abn:::buildScoreCache(data.df=df, data.dists=dist, max.parents=3,method="mle")
    mycache.old.1.1 <- abn:::buildScoreCache(data.df=df, data.dists=dist, max.parents=3, dry.run=TRUE, method="bayes")
    ## test wrapper buildScoreCache method
    expect_equal(mycache.computed.mle.1,unclass(mycache.computed.mle.1.1))
    expect_equal(mycache.old.1,unclass(mycache.old.1.1))
    ## test

    ## dag retain
    mycache.computed.mle.2 <- abn:::buildScoreCache.mle(data.df=df, data.dists=dist, dag.banned=NULL, dag.retained=~Ozone |
        Solar.R, max.parents=3, dry.run=TRUE)
 #   if(requireNamespace("INLA", quietly=TRUE)){
    mycache.old.2 <- abn:::buildScoreCache.bayes(data.df=df, data.dists=dist, max.parents=3, dry.run=TRUE, dag.retained=~Ozone |
        Solar.R, dag.banned=NULL)
 #   }  else expect_equal(1,1)

    mycache.computed.mle.3 <- abn:::buildScoreCache.mle(data.df=df, data.dists=dist, dag.banned=NULL, dag.retained=~Wind |
        ., max.parents=6, dry.run=TRUE)
    mycache.old.3 <- abn:::buildScoreCache.bayes(data.df=df, data.dists=dist, max.parents=6, dry.run=TRUE, dag.retained=~Wind |
        ., dag.banned=NULL)


    ## dag ban
    mycache.computed.mle.4 <- abn:::buildScoreCache.mle(data.df=df, data.dists=dist, dag.banned=~Ozone | Solar.R:Wind,
        max.parents=3, dry.run=TRUE)
    mycache.old.4 <- abn:::buildScoreCache.bayes(data.df=df, data.dists=dist, max.parents=3, dry.run=TRUE, dag.banned=~Ozone |
        Solar.R:Wind)

    ## test cache
 #   if(requireNamespace("INLA", quietly=TRUE)){
    expect_that(mycache.computed.mle$children, equals(mycache.old$children))
    expect_that(mycache.computed.mle$node.defn, equals(mycache.old$node.defn))
    expect_that(mycache.computed.mle.1$children, equals(mycache.old.1$children))
    expect_that(mycache.computed.mle.1$node.defn, equals(mycache.old.1$node.defn))
    expect_that(mycache.computed.mle.2$children, equals(mycache.old.2$children))
    expect_that(mycache.computed.mle.2$node.defn, equals(mycache.old.2$node.defn))
 #   } else expect_equal(1,1)

    expect_that(mycache.computed.mle.3$children, equals(mycache.old.3$children))
    expect_that(mycache.computed.mle.3$node.defn, equals(mycache.old.3$node.defn))
    expect_that(mycache.computed.mle.4$children, equals(mycache.old.4$children))
    expect_that(mycache.computed.mle.4$node.defn, equals(mycache.old.4$node.defn))

    ## which.nodes
    expect_error(abn:::buildScoreCache.mle(data.df=df, data.dists=dist[1:5], max.parents=6, which.nodes=1:5), NA)

    ## max.parents
    expect_equal(abn:::buildScoreCache.mle(data.df=df, data.dists=dist[1:5], max.parents=6, which.nodes=1:5), abn:::buildScoreCache.mle(data.df=df,
        data.dists=dist[1:5], max.parents=list(Ozone=4, Solar.R=4, Wind=4, Temp=4, Month=4), which.node=1:5))

    ## test scoring Gaussian
    dist <- list(a="gaussian", b="gaussian")

    data.param <- matrix(data=c(0, 0.5, 0, 0), nrow=2L, ncol=2L, byrow=TRUE)

    # naming
    colnames(data.param) <- rownames(data.param) <- names(dist)

    out.sim <- invisible(simulateAbn(data.dists=dist, n.chains=1, n.adapt=100, n.thin=1, n.iter=100, data.param=data.param,
        simulate=TRUE, seed=132,verbose=FALSE))

    mycache <- invisible(abn:::buildScoreCache.mle(data.df=out.sim, data.dists=dist, max.parents=2, centre=FALSE))

    # mLik
    expect_that(mycache$mlik[1], equals(as.numeric(logLik(lm(formula=out.sim$a ~ 1)))))
    expect_that(mycache$mlik[2], equals(as.numeric(logLik(lm(formula=out.sim$a ~ 1 + out.sim$b)))))
    expect_that(mycache$mlik[3], equals(as.numeric(logLik(lm(formula=out.sim$b ~ 1)))))

    # AIC
    expect_that(mycache$aic[1], equals(as.numeric(AIC(lm(formula=out.sim$a ~ 1)))))
    expect_that(mycache$aic[2], equals(as.numeric(AIC(lm(formula=out.sim$a ~ 1 + out.sim$b)))))
    expect_that(mycache$aic[3], equals(as.numeric(AIC(lm(formula=out.sim$b ~ 1)))))

    # BIC
    expect_that(mycache$bic[1], equals(as.numeric(BIC(lm(formula=out.sim$a ~ 1)))))
    expect_that(mycache$bic[2], equals(as.numeric(BIC(lm(formula=out.sim$a ~ 1 + out.sim$b)))))
    expect_that(mycache$bic[3], equals(as.numeric(BIC(lm(formula=out.sim$b ~ 1)))))

    ## Gaussian: correlation coefficients

    mycache <- invisible(abn:::buildScoreCache.mle(data.df=out.sim, data.dists=dist, max.parents=2, centre=TRUE))

    expect_that(mycache$mlik[1], equals(as.numeric(logLik(lm(formula=(out.sim$a - mean(out.sim$a))/sd(out.sim$a) ~ 1)))))

    # Binomial
    dist <- list(a="binomial", b="binomial")

    data.param <- matrix(data=c(0, 0.5, 0, 0), nrow=2L, ncol=2L, byrow=TRUE)

    # naming
    colnames(data.param) <- rownames(data.param) <- names(dist)

    out.sim <- invisible(simulateAbn(data.dists=dist, n.chains=1, n.adapt=100, n.thin=1, n.iter=100, data.param=data.param,
        simulate=TRUE, seed=132,verbose=FALSE))

    mycache <- invisible(abn:::buildScoreCache.mle(data.df=out.sim, data.dists=dist, max.parents=2, centre=FALSE))

    # mLik
    expect_that(mycache$mlik[1], equals(as.numeric(logLik(glm(formula=out.sim$a ~ 1, family=binomial)))))
    expect_that(mycache$mlik[2], equals(as.numeric(logLik(glm(formula=out.sim$a ~ 1 + out.sim$b, family=binomial)))))
    expect_that(mycache$mlik[3], equals(as.numeric(logLik(glm(formula=out.sim$b ~ 1, family=binomial)))))

    # AIC
    expect_that(mycache$aic[1], equals(as.numeric(AIC(glm(formula=out.sim$a ~ 1, family=binomial)))))
    expect_that(mycache$aic[2], equals(as.numeric(AIC(glm(formula=out.sim$a ~ 1 + out.sim$b, family=binomial)))))
    expect_that(mycache$aic[3], equals(as.numeric(AIC(glm(formula=out.sim$b ~ 1, family=binomial)))))

    # BIC
    expect_that(mycache$bic[1], equals(as.numeric(BIC(glm(formula=out.sim$a ~ 1, family=binomial)))))
    expect_that(mycache$bic[2], equals(as.numeric(BIC(glm(formula=out.sim$a ~ 1 + out.sim$b, family=binomial)))))
    expect_that(mycache$bic[3], equals(as.numeric(BIC(glm(formula=out.sim$b ~ 1, family=binomial)))))

    # Poisson
    dist <- list(a="poisson", b="poisson")

    data.param <- matrix(data=c(0, 0.5, 0, 0), nrow=2L, ncol=2L, byrow=TRUE)

    # naming
    colnames(data.param) <- rownames(data.param) <- names(dist)

    out.sim <- invisible(simulateAbn(data.dists=dist, n.chains=1, n.adapt=100, n.thin=1, n.iter=100, data.param=data.param,
        simulate=TRUE, seed=132,verbose=FALSE))

    mycache <- invisible(abn:::buildScoreCache.mle(data.df=out.sim, data.dists=dist, max.parents=2))

    # mLik
    expect_that(mycache$mlik[1], equals(as.numeric(logLik(glm(formula=out.sim$a ~ 1, family=poisson)))))
    expect_that(mycache$mlik[2], equals(as.numeric(logLik(glm(formula=out.sim$a ~ 1 + out.sim$b, family=poisson)))))
    expect_that(mycache$mlik[3], equals(as.numeric(logLik(glm(formula=out.sim$b ~ 1, family=poisson)))))

    # AIC
    expect_that(mycache$aic[1], equals(as.numeric(AIC(glm(formula=out.sim$a ~ 1, family=poisson)))))
    expect_that(mycache$aic[2], equals(as.numeric(AIC(glm(formula=out.sim$a ~ 1 + out.sim$b, family=poisson)))))
    expect_that(mycache$aic[3], equals(as.numeric(AIC(glm(formula=out.sim$b ~ 1, family=poisson)))))

    # BIC
    expect_that(mycache$bic[1], equals(as.numeric(BIC(glm(formula=out.sim$a ~ 1, family=poisson)))))
    expect_that(mycache$bic[2], equals(as.numeric(BIC(glm(formula=out.sim$a ~ 1 + out.sim$b, family=poisson)))))
    expect_that(mycache$bic[3], equals(as.numeric(BIC(glm(formula=out.sim$b ~ 1, family=poisson)))))

    # multinomial (as response)
    dist <- list(a="multinomial", b="gaussian")

    # link matrix
    data.param <- matrix(data=c(1, 0, 0, 1), nrow=2L, ncol=2L, byrow=TRUE)
    # naming
    colnames(data.param) <- rownames(data.param) <- names(dist)

    data.param.mult <- matrix(data=c(0.7, 0.1, 0.2, 0, 0, 0), nrow=2L, ncol=3L, byrow=TRUE)

    out <- invisible(simulateAbn(data.dists=dist, n.chains=1, n.adapt=100, n.thin=1, n.iter=1000, data.param=data.param,
        simulate=TRUE, data.param.mult=data.param.mult, seed=132,verbose=FALSE))

    res <- out$`a[1]`
    res[out$`a[2]` == 1] <- 2
    res[out$`a[3]` == 1] <- 3

    dta <- data.frame(b=rnorm(n=length(res), mean=1, sd=1), c=rnorm(n=length(res), mean=5, sd=1))

    dta$res <- factor(res)

    dist <- list(a="gaussian", b="gaussian", c="multinomial")
suppressWarnings(require(boot))

    names(dta) <- names(dist)

    mycache <- invisible(abn:::buildScoreCache.mle(data.df=dta, data.dists=dist, max.parents=3, centre=FALSE))

    # mLik
    expect_that(mycache$mlik[9], equals(as.numeric(logLik(multinom(formula=dta$c ~ 1, trace=FALSE)))))
    expect_that(mycache$mlik[10], equals(as.numeric(logLik(multinom(formula=dta$c ~ 1 + dta$a, trace=FALSE)))))
    expect_that(mycache$mlik[11], equals(as.numeric(logLik(multinom(formula=dta$c ~ 1 + dta$b, trace=FALSE)))))

    # AIC
    expect_that(mycache$aic[9], equals(as.numeric(AIC(multinom(formula=dta$c ~ 1, trace=FALSE)))))
    expect_that(mycache$aic[10], equals(as.numeric(AIC(multinom(formula=dta$c ~ 1 + dta$a, trace=FALSE)))))
    expect_that(mycache$aic[11], equals(as.numeric(AIC(multinom(formula=dta$c ~ 1 + dta$b, trace=FALSE)))))

    # BIC
    expect_that(mycache$bic[9], equals(as.numeric(BIC(multinom(formula=dta$c ~ 1, trace=FALSE)))))
    expect_that(mycache$bic[10], equals(as.numeric(BIC(multinom(formula=dta$c ~ 1 + dta$a, trace=FALSE)))))
    expect_that(mycache$bic[11], equals(as.numeric(BIC(multinom(formula=dta$c ~ 1 + dta$b, trace=FALSE)))))

    ## Multinomial (as predictors)

    # mLik
    m.1 <- model.matrix(object=~dta$c + 0)#, contrasts.arg=contr.treatment)

    expect_that(mycache$mlik[3], equals(as.numeric(logLik(glm(formula=dta$a ~ 0 + m.1, family=gaussian)))))
    expect_that(mycache$mlik[4], equals(as.numeric(logLik(glm(formula=dta$a ~ 0 + m.1 + dta$b, family=gaussian)))))

    # AIC
    expect_that(mycache$aic[3], equals(as.numeric(AIC(glm(formula=dta$a ~ 0 + m.1, family=gaussian)))))
    expect_that(mycache$aic[4], equals(as.numeric(AIC(glm(formula=dta$a ~ 0 + m.1 + dta$b, family=gaussian)))))

    # BIC
    expect_that(mycache$bic[3], equals(as.numeric(BIC(glm(formula=dta$a ~ 0 + m.1, family=gaussian)))))
    expect_that(mycache$bic[4], equals(as.numeric(BIC(glm(formula=dta$a ~ 0 + m.1 + dta$b, family=gaussian)))))

    # data separation

    ## simulation data
    n <- 1000
    x1 <- rnorm(n)
    x2 <- rbinom(n, 1, 0.5)
    b0 <- 1
    b1 <- 1.5
    b2 <- 2
    y <- rbinom(n, 1, boot::inv.logit(b0 + b1 * x1 + b2 * x2))

    y <- ifelse(x2 == 1, 1, y)

    dist <- list(a="binomial", b="gaussian", c="binomial")
    dta <- data.frame(y, x1, x2)
    names(dta) <- names(dist)

# suppressWarnings(require(brglm))

#    mycache <- invisible(abn:::buildScoreCache.mle(data.df=dta, data.dists=dist, max.parents=2, centre=FALSE, dry.run=TRUE))  ##
#   save(mycache, file='testdata/mycache.Rdata')
     if(requireNamespace("brglm", quietly=TRUE)){
       load(file='testdata/mycache.Rdata')

      expect_that(mycache$mlik[1],  equals(suppressWarnings(as.numeric(logLik(brglm::brglm(formula=dta$a ~ 1)))), tolerance=0.01))
      expect_that(mycache$mlik[2],  equals(as.numeric(logLik(suppressWarnings(brglm::brglm(formula=dta$a ~ dta$b)))), tolerance=0.01))
      expect_that(mycache$mlik[3],  equals(as.numeric(logLik(suppressWarnings(brglm::brglm(formula=dta$a ~ dta$c)))), tolerance=0.01))
      expect_that(mycache$mlik[4],  equals(as.numeric(logLik(suppressWarnings(brglm::brglm(formula=dta$a ~ dta$b + dta$c)))), tolerance=0.01))
      expect_that(mycache$mlik[12], equals(as.numeric(logLik(suppressWarnings(brglm::brglm(formula=dta$c ~ dta$b + dta$a)))), tolerance=0.01))
     }  else {
        cat('Package brglm not available, number of passed tests might be different')
    }
# else expect_equal(1,1)
})

test_that("Markov Blanket", {
    dist <- list(a="gaussian", b="gaussian", c="gaussian", d="gaussian", e="binomial", f="binomial")

    # define parameter matrix
    data.param <- matrix(data=c(0, 0.2, 0.5, 0, 0.01, 0, 0, 0, 0.3, 0.1, 0, 0.8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0.1, 0, 0, 0, 0, 0, 0), nrow=6L, ncol=6L, byrow=TRUE)
    colnames(data.param) <- rownames(data.param) <- names(dist)

    a <- mb(dag=data.param, node="b", data.dists=dist)
    b <- mb(dag=data.param, node="e", data.dists=dist)
    c <- mb(dag=data.param, node=c("b", "e"), data.dists=dist)

    expect_that(a, equals(c("a", "c", "d", "f", "e")))
    expect_that(b, equals(c("a", "f", "b", "c")))
    expect_that(c, equals(c("a", "c", "d", "f", "e", "b")))
})

test_that("Formula statement", {
    dist <- list(a="gaussian", b="gaussian", c="gaussian", d="gaussian", e="gaussian", f="gaussian")

    m.formula.1 <- createAbnDag(dag=~a | b:c + b | c:d + a | e:f, data.dists=dist)$dag
    m.formula.2 <- createAbnDag(dag=~a | ., data.dists=dist)$dag

    m.true.1 <- matrix(data=c(0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0), nrow=6, ncol=6, byrow=TRUE)
    m.true.2 <- matrix(data=c(0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0), nrow=6, ncol=6, byrow=TRUE)
    colnames(m.true.1) <- rownames(m.true.1) <- colnames(m.true.2) <- rownames(m.true.2) <- names(dist)

    expect_that(m.formula.1, equals(m.true.1))
    expect_that(m.formula.2, equals(m.true.2))

    ## formula with real data

    df <- airquality[complete.cases(airquality), ]

    # distribution (gaussian)
    dist <- list(Ozone="gaussian", Solar.R="gaussian", Wind="gaussian", Temp="gaussian", Month="gaussian", Day="gaussian")
    names(dist) <- colnames(df)

    m.formula.1 <- createAbnDag(dag=~Ozone | Solar.R, data.dists=dist)$dag
    m.formula.2 <- createAbnDag(dag=~Solar.R | ., data.dists=dist)$dag

    m.true.1 <- matrix(data=c(0, 1, 0, 0, 0, 0, rep(0, 30)), nrow=6, ncol=6, byrow=TRUE)
    m.true.2 <- matrix(data=c(0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, rep(0, 24)), nrow=6, ncol=6, byrow=TRUE)
    colnames(m.true.1) <- rownames(m.true.1) <- colnames(m.true.2) <- rownames(m.true.2) <- names(dist)

    expect_that(m.formula.1, equals(m.true.1))
    expect_that(m.formula.2, equals(m.true.2))


})

test_that("infoDag()", {

    dag <- matrix(data=0, nrow=6, ncol=6)
    dist <- list(a="gaussian", b="gaussian", c="gaussian", d="gaussian", e="gaussian", f="gaussian")
    colnames(dag) <- rownames(dag) <- names(dist)

    infoDag.out1 <- infoDag(dag=dag,node.names=names(dist))

    dag <- matrix(data=c(0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0), nrow=6, ncol=6)
    colnames(dag) <- rownames(dag) <- names(dist)

    infoDag.out2 <- infoDag(dag=dag, node.names=names(dist))

    dag <- matrix(data=c(0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0), nrow=6, ncol=6)
    colnames(dag) <- rownames(dag) <- names(dist)

    infoDag.out3 <- infoDag(dag=dag, node.names=names(dist))

    expect_equal(infoDag.out1$n.nodes, 6)
    expect_equal(infoDag.out1$n.arcs, 0)
    expect_equal(infoDag.out1$mb.average, 0)
    expect_equal(infoDag.out1$nh.average, 0)
    expect_equal(infoDag.out1$parent.average, 0)
    expect_equal(infoDag.out1$children.average, 0)

    expect_equal(infoDag.out2$n.nodes, 6)
    expect_equal(infoDag.out2$n.arcs, 2)
    expect_equal(infoDag.out2$mb.average, 1)
    expect_equal(infoDag.out2$nh.average, 2/3)
    expect_equal(infoDag.out2$parent.average, 1/3)
    expect_equal(infoDag.out2$children.average, 1/3)

    expect_equal(infoDag.out3$n.nodes, 6)
    expect_equal(infoDag.out3$n.arcs, 3)
    expect_equal(infoDag.out3$mb.average, 2)
    expect_equal(infoDag.out3$nh.average, 1)
    expect_equal(infoDag.out3$parent.average, 0.5)
    expect_equal(infoDag.out3$children.average, 0.5)

})

test_that("logit(), expit()", {

    expect_equal(abn::logit(x=0.678), boot::logit(p=0.678))
    expect_equal(abn::logit(x=0.783491741), boot::logit(p=0.783491741))

    expect_equal(abn::expit(x=0.678), boot::inv.logit(x=0.678))
    expect_equal(abn::expit(x=-0.783492343421741), boot::inv.logit(x=-0.783492343421741))

})

test_that("discretization(), entropyData()", {

    suppressWarnings(require(entropy))

    dist <- list(a="gaussian", b="gaussian", c="gaussian")
    data.param <- matrix(data=c(0, 1, 0, 0, 0, 1, 0, 0, 0), nrow=3L, ncol=3L, byrow=TRUE)

    data.param.var <- matrix(data=0, nrow=3L, ncol=3L)
    diag(data.param.var) <- c(0.1, 0.1, 0.1)

    out <- invisible(simulateAbn(data.dists=dist, n.chains=1, n.adapt=1e+06, n.thin=1, n.iter=100, data.param=data.param,
        data.param.var=data.param.var, simulate=TRUE, seed=132,verbose=FALSE))

    y2d.entropy.1 <- entropy::discretize2d(x1=out$a, x2=out$b, numBins1=100, numBins2=100)
    y2d.abn.1 <- abn::discretization(data.df=out[, c(1, 2)], discretization.method=100, data.dists=dist[c(1, 2)], nb.states=FALSE)

    y2d.entropy.2 <- entropy::discretize2d(x1=out$a, x2=out$c, numBins1=100, numBins2=100)
    y2d.abn.2 <- abn::discretization(data.df=out[, c(1, 3)], discretization.method=100, data.dists=dist[c(1, 2)], nb.states=FALSE)


    # Not the same dimnames!
    dimnames(y2d.abn.1) <- dimnames(y2d.entropy.1)
    dimnames(y2d.abn.2) <- dimnames(y2d.entropy.2)

    expect_equal(y2d.abn.1, y2d.entropy.1)
    expect_equal(y2d.abn.2, y2d.entropy.2)

    expect_equal(entropyData(freqs.table=y2d.abn.1), entropy.empirical(y=y2d.entropy.1, unit="log2"))
    expect_equal(entropyData(freqs.table=y2d.abn.2), entropy.empirical(y=y2d.entropy.2, unit="log2"))
})

test_that("or(), miData(), ", {

    prob <- logit(0.75)

    ## pure binomial
    dist <- list(a="binomial", b="binomial", c="binomial")
    data.param <- matrix(data=c(0, prob, 0, 0, 0, prob, 0, 0, 0), nrow=3L, ncol=3L, byrow=TRUE)

    out <- invisible(simulateAbn(data.dists=dist, n.chains=1, n.adapt=1e+06, n.thin=10, n.iter=100, data.param=data.param,
        simulate=TRUE, seed=132,verbose=FALSE))

    x.bc <- table(as.numeric(as.character(out$b)), as.numeric(as.character(out$c)), dnn=c("b", "c"))

    x.ab <- table(as.numeric(out$a), as.numeric(out$b), dnn=c("a", "b"))

    x.ac <- table(as.numeric(out$a), as.numeric(out$c), dnn=c("a", "c"))

    ## OR()

    expect_equal(abn::or(x=x.bc), 3/2)
    expect_equal(abn::or(x=x.ab), 1)
    expect_equal(abn::or(x=x.ac), 8/3)

    ## MUTUAL INFORMATION

    expect_equal(miData(freqs.table=x.bc, method="mi.raw"), mi.empirical(y2d=x.bc))
    expect_equal(miData(freqs.table=x.ab, method="mi.raw"), mi.empirical(y2d=x.ab))
    expect_equal(miData(freqs.table=x.ac, method="mi.raw"), mi.empirical(y2d=x.ac))

})

test_that("skewness()", {

    suppressWarnings(require(moments))

    data <- c(19.09, 19.55, 17.89, 17.73, 25.15, 27.27, 25.24, 21.05, 21.65, 20.92, 22.61, 15.71, 22.04, 22.6, 24.25)

    expect_equal(abn:::skewness(x=data), moments::skewness(x=data))

})

test_that("essentialGraph()", {

    dist.test <- list(a="gaussian", b="gaussian", c="gaussian")

    ## essentialGraph()
    expect_equal(object=essentialGraph(  dag=~a | b + c | b, node.names=names(dist.test)),
                 expected=essentialGraph(dag=~b | a + c | b, node.names=names(dist.test)))

    expect_equal(object=  essentialGraph(dag=~a | b + b | c, node.names=names(dist.test)),
                 expected=essentialGraph(dag=~b | a + c | b, node.names=names(dist.test)))

    # more complex
    dist.test <- list(a="gaussian", b="gaussian", c="gaussian", d="gaussian", e="gaussian", f="gaussian")
    # examples from 'Bayesian Networks in Bioinformatics, Kyu-Baek Hwang'

    minimal.dag <- matrix(data=c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
        1, 0, 0, 0, 1, 0), nrow=6, byrow=TRUE)
    colnames(minimal.dag) <- rownames(minimal.dag) <- names(dist.test)
    completed.dag <- matrix(data=c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
        1, 0, 0, 0, 1, 0), nrow=6, byrow=TRUE)
    colnames(completed.dag) <- rownames(completed.dag) <- names(dist.test)

    expect_equal(essentialGraph(dag=~f | e:a + e | c + c | a:b + d | b, node.names=names(dist.test),
        PDAG="minimal"), minimal.dag)
    expect_equal(essentialGraph(dag=~f | e:a + e | c + c | a:b + d | b, node.names=names(dist.test),
        PDAG="completed"),  completed.dag)

})

test_that("compareDag()", {

    a <- matrix(data=0, nrow=3, ncol=3)

    a1 <- matrix(data=c(0, 0, 0, 1, 0, 0, 1, 0, 0), nrow=3, ncol=3)
    a2 <- matrix(data=c(0, 0, 0, 1, 0, 0, 1, 1, 0), nrow=3, ncol=3)
    b <- matrix(data=c(0, 0, 0, 1, 0, 1, 1, 0, 0), nrow=3, ncol=3)

    expect_equal(suppressWarnings(compareDag(ref=a, test=b))$`Hamming-distance`, expected=3)
    expect_equal(compareDag(ref=a1, test=b)$`Hamming-distance`, expected=1)
    expect_equal(compareDag(ref=a1, test=b)$TPR, expected=1)
    expect_equal(compareDag(ref=a1, test=b)$PPV, expected=2/3)
    expect_equal(compareDag(ref=a2, test=b)$`Hamming-distance`, expected=1)
    expect_equal(compareDag(ref=a2, test=b)$PPV, expected=2/3)
    expect_equal(compareDag(ref=a2, test=b)$FDR, expected=1/3)
    expect_equal(compareDag(ref=a2, test=b)$TPR, expected=2/3)
})

test_that("scoreContribution()", {

    mydat <- ex1.dag.data[,c("b1","g1","p1")]
    ## take a subset of cols

    ## setup distribution list for each node
    mydists <- list(b1="binomial",
                    g1="gaussian",
                    p1="poisson"
    )

    ## now build cache
    mycache <- buildScoreCache(data.df=mydat,data.dists=mydists,max.parents=1, method="mle")

    #now find the globally best DAG
    mp.dag <- mostProbable(score.cache=mycache, score="bic", verbose=FALSE)

    out <- scoreContribution(object=mp.dag)

    out.fit <- fitAbn(object=mp.dag, method="mle")

    expect_equal(unname(colSums(out$mlik))/1000, unname(unlist(out.fit$mliknode))/1000,tolerance=1e-2)

})

test_that("plotabn()", {

    #Define distribution list
    dist <- list(a="gaussian", b="gaussian", c="gaussian", d="gaussian", e="binomial", f="multinomial")

    #Plot from a formula and markov blanket for multinomial node
    expect_warning(plotabn(dag.m=~a|b:c:e+b|c:d:f+e|f,markov.blanket.node="b", data.dists =dist))

    dist <- list(a="gaussian", b="gaussian")

    data.param <- matrix(data=c(0, 0.5, 0, 0), nrow=2L, ncol=2L, byrow=TRUE)

    # naming
    colnames(data.param) <- rownames(data.param) <- names(dist)

    out.sim <- invisible(simulateAbn(data.dists=dist, n.chains=1, n.adapt=100, n.thin=1, n.iter=100, data.param=data.param,
                                     simulate=TRUE, seed=132,verbose=FALSE))

    mycache <- invisible(abn:::buildScoreCache.mle(data.df=out.sim, data.dists=dist, max.parents=2, centre=FALSE))
    mycache <- buildScoreCache(data.df=out.sim, data.dists=dist, max.parents=2, centre=FALSE, method='mle')

    dag <- mostProbable(score.cache=mycache, verbose=FALSE)

    # class abnlearned
    expect_silent(plotAbn(dag=dag))

})
