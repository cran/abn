###############################################################################
## fitabn.R ---
## Author : Fraser Lewis & Gilles Kratzer
## Last modified : 02/12/2012
## Last modified : 29/09/2014 by Marta Pittavino, just renamed.
##               :  06/12/2016 by GK (implementation of formula
############################################################################### statment)

## fit a given DAG to data

fitabn <- function(...) {
    .Deprecated("fitAbn", msg="'fitabn' is deprecated.\n Use 'fitAbn' instead but note that arguments have slightly changed.")
    fitAbn(...)
}

## control function

fit.control <- function (method = "bayes", mean = 0, prec = 0.001, loggam.shape = 1,
                         loggam.inv.scale = 5e-05, max.mode.error = 10, max.iters = 100,
                         epsabs = 1e-07, error.verbose = FALSE, trace=0L, epsabs.inner = 1e-06,
                         max.iters.inner = 100, finite.step.size = 1e-07, hessian.params = c(1e-04, 0.01),
                         max.iters.hessian = 10, max.hessian.error = 1e-04, factor.brent = 100,
                         maxiters.hessian.brent = 10, num.intervals.brent = 100, min.pdf = 0.001,
                         n.grid = 250, std.area = TRUE, marginal.quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975),
                         max.grid.iter = 1000, marginal.node = NULL, marginal.param = NULL, variate.vec = NULL,
                         max.irls = 100, tol = 10^-11, seed = 9062019)
  {

  if(method == "bayes"){

    ctrl <- list(mean = mean, prec = prec, loggam.shape = loggam.shape, loggam.inv.scale = loggam.inv.scale,
         max.mode.error = max.mode.error, max.iters = max.iters, epsabs = epsabs,
         error.verbose = error.verbose, trace=trace, epsabs.inner = epsabs.inner, max.iters.inner = max.iters.inner,
         finite.step.size = finite.step.size, hessian.params = hessian.params,
         max.iters.hessian = max.iters.hessian, max.hessian.error = max.hessian.error,
         factor.brent = factor.brent, maxiters.hessian.brent = maxiters.hessian.brent,
         num.intervals.brent = num.intervals.brent, min.pdf = min.pdf, n.grid = n.grid,
         std.area = std.area, marginal.quantiles = marginal.quantiles,
         max.grid.iter = max.grid.iter, marginal.node = marginal.node, marginal.param = marginal.param,
         variate.vec = variate.vec, seed = seed)
  }

  if(method == "mle"){

    ctrl <- list(max.irls = max.irls, tol = tol, seed = seed)
  }

  return(ctrl)

}


## fit a given DAG to data
fitAbn <- function(object = NULL, dag = NULL, data.df = NULL, data.dists = NULL,
    method = NULL, group.var = NULL, adj.vars = NULL, cor.vars = NULL, centre = TRUE,
    compute.fixed = FALSE, control=NULL, verbose = FALSE, ...) {

    if(inherits(x = dag, what = "abnLearned") & is.null(object)){
        object <- dag
        message("Best practice with abn > 2.0 requires to pass 'dag' as 'object' parameter.")
    }

    ## method abnCache
    if (!is.null(object) & inherits(x=object, what="abnLearned")) {

        dag <- object$dag
        object <- object$score.cache

        data.df <- object$data.df
        data.dists <- object$data.dists
        # group.var <- object$group.var
        # cor.vars <- object$cor.vars
        # adj.vars <- object$adj.vars

        fitmethod <- object$method
    } else  fitmethod <- NULL

    if (is.null(method))  {
        method <- if (is.null(fitmethod)) "bayes" else fitmethod
    }else{
        if (verbose & (!is.null(fitmethod)))
            if (method != fitmethod) cat("Fitting and learned methods differ!\n")
    }


    ## start tests
    method <- c("bayes", "mle")[pmatch(tolower(method), c("bayes", "mle"))][1]
    if (is.na(method)) stop("'method' should be 'bayes' or 'mle'.")

    ctrl <- fit.control(method = method)
    if (!missing(control)) {
      control <- as.list(control)
      ctrl[names(control)] <- control
    }

    nmsC <- names(ctrl)
    ctrl[(namc <- names(control))] <- control
    if (length(noNms <- namc[!namc %in% nmsC]))
        warning("Unknown names in 'control': ", paste(noNms, collapse=", "))


    if (method == "bayes") {

        out <- fitAbn.bayes(dag, data.df = data.df, data.dists = data.dists,
                            group.var = group.var, cor.vars = cor.vars,
                            compute.fixed = compute.fixed, control = ctrl)

        return(out)
    }
    if (method == "mle") {

        out <- fitAbn.mle(dag, data.df = data.df, data.dists = data.dists, adj.vars = adj.vars,
                cor.vars = cor.vars, verbose = verbose, centre = centre, control = ctrl)

        return(out)
    }

}  #EOF
