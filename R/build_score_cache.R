############################################################################### buildscorecache.R

buildscorecache <- function(...) {
    .Deprecated("buildScoreCache")#, msg="'mostprobable' is deprecated.\n Use 'mostProbable' instead but note that arguments have slightly changed.")
    buildScoreCache(...)
}

build.control <- function (method = "bayes", max.mode.error = 10, mean = 0, prec = 0.001,
                           loggam.shape = 1, loggam.inv.scale = 5e-05, max.iters = 100,
                           epsabs = 1e-07, error.verbose = FALSE, trace = 0L, epsabs.inner = 1e-06,
                           max.iters.inner = 100, finite.step.size = 1e-07, hessian.params = c(1e-04, 0.01),
                           max.iters.hessian = 10, max.hessian.error = 0.5, factor.brent = 100,
                           maxiters.hessian.brent = 100, num.intervals.brent = 100,
                           ncores = 0, max.irls = 100, tol = 10^-8, seed = 9062019) {

  if(method == "bayes"){
    ctrl <- list(max.mode.error = max.mode.error, mean = mean, prec = prec, loggam.shape = loggam.shape,
         loggam.inv.scale = loggam.inv.scale, max.iters = max.iters, epsabs = epsabs,
         error.verbose = error.verbose, trace = trace, epsabs.inner = epsabs.inner, max.iters.inner = max.iters.inner,
         finite.step.size = finite.step.size, hessian.params = hessian.params,
         max.iters.hessian = max.iters.hessian, max.hessian.error = max.hessian.error,
         factor.brent = factor.brent, maxiters.hessian.brent = maxiters.hessian.brent,
         num.intervals.brent = num.intervals.brent, seed = seed)
  }

  if(method == "mle"){
    ctrl <- list(ncores = 0, max.iters = max.iters, tol = tol, seed = seed)
  }
  return(ctrl)
}

## precompute a cache of scores to data
buildScoreCache <- function(data.df = NULL, data.dists = NULL, method = "bayes",
                            group.var = NULL, adj.vars = NULL, cor.vars = NULL,
                            dag.banned = NULL, dag.retained = NULL, max.parents = NULL,
                            which.nodes = NULL, defn.res = NULL, centre = TRUE,
                            dry.run = FALSE, control = NULL, verbose = FALSE,
                            ...) {

    ## start tests
    method <- tolower(method)
    method <- c("bayes", "mle")[pmatch(method, c("bayes", "mle"))]
    if (is.na(method)) stop("'method' should be 'bayes' or 'mle'.")

    data.dists <- validate_dists( data.dists, returnDists=TRUE)

    ctrl <- build.control(method = method)
    if (!missing(control)) {
      control <- as.list(control)
      ctrl[names(control)] <- control
    }

    nmsC <- names(ctrl)

    ctrl[(namc <- names(control))] <- control
    if (length(noNms <- namc[!namc %in% nmsC]))
        warning("Unknown names in 'control': ", paste(noNms, collapse = ", "))



    if (method == "bayes") {

        out <- buildScoreCache.bayes(data.df = data.df, data.dists = data.dists,
                                     group.var = group.var, cor.vars = cor.vars,
                                     dag.banned = dag.banned, dag.retained = dag.retained,
                                     max.parents = max.parents, which.nodes = which.nodes,
                                     defn.res = defn.res, dry.run = dry.run,
                                     verbose = verbose, centre = centre,
                                     control = ctrl)

    }
    if (method == "mle") {

        if (!is.null(group.var)) warning("'group.var' is ignored with method 'mle'.")

            
        out <- buildScoreCache.mle(data.df = data.df, data.dists = data.dists,
                                   adj.vars = adj.vars, cor.vars = cor.vars,
                                   dag.banned = dag.banned, dag.retained = dag.retained,
                                   max.parents = max.parents,
                                   which.nodes = which.nodes, defn.res = defn.res,
                                   dry.run = dry.run, verbose = verbose, centre = centre,
                                   control = ctrl)

    }
    class(out) <- c("abnCache")
    return(out)
    
}  #EOF




